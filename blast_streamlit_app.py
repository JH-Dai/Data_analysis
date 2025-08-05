import streamlit as st
import pandas as pd
import os
import io

st.set_page_config(page_title="BLAST筛选工具", layout="wide")
st.title("🧬 BLAST 结果分析工具（上传 + 拆分 + 筛选 + 下载）")

UPLOAD_DIR = "uploaded"
SPLIT_DIR = "split_queries"
os.makedirs(UPLOAD_DIR, exist_ok=True)
os.makedirs(SPLIT_DIR, exist_ok=True)

# 会话状态
if "results_ready" not in st.session_state:
    st.session_state.results_ready = False
if "global_result_df" not in st.session_state:
    st.session_state.global_result_df = pd.DataFrame()
if "total_hits" not in st.session_state:
    st.session_state.total_hits = 0
if "skipped_queries" not in st.session_state:
    st.session_state.skipped_queries = []

# 上传 + 拆分
uploaded_file = st.file_uploader("📤 请上传 BLASTN 结果文件（含 `# Query:` 行）", type=["txt", "tsv"])

def split_file(file_path):
    with open(file_path, "r", encoding="utf-8") as f:
        current_block = []
        query_id = None
        block_count = 0
        for line in f:
            if line.startswith("# Query:"):
                if current_block and query_id:
                    out_path = os.path.join(SPLIT_DIR, f"{query_id}.txt")
                    with open(out_path, "w", encoding="utf-8") as out:
                        out.writelines(current_block)
                    block_count += 1
                query_id = line.strip().split(":")[1].strip().replace("/", "_")
                current_block = [line]
            else:
                current_block.append(line)
        if current_block and query_id:
            out_path = os.path.join(SPLIT_DIR, f"{query_id}.txt")
            with open(out_path, "w", encoding="utf-8") as out:
                out.writelines(current_block)
            block_count += 1
    return block_count

if uploaded_file is not None:
    file_path = os.path.join(UPLOAD_DIR, uploaded_file.name)
    with open(file_path, "wb") as f:
        f.write(uploaded_file.getbuffer())
    st.success(f"✅ 文件已上传：{uploaded_file.name}")
    count = split_file(file_path)
    st.success(f"🔧 已成功拆分为 {count} 个 query blocks，存放于 `{SPLIT_DIR}/`")

# 查询列表
query_files = sorted([f for f in os.listdir(SPLIT_DIR) if f.endswith(".txt")])
if not query_files:
    st.info("📭 当前无 query 数据，请先上传 BLAST 结果文件。")
    st.stop()

# 读取 query 文件
def load_query_file(filepath):
    with open(filepath, 'r', encoding="utf-8") as f:
        lines = f.readlines()
    data_lines = [line for line in lines if not line.startswith("#") and line.strip()]
    if not data_lines:
        return pd.DataFrame(), lines
    df = pd.read_csv(io.StringIO("".join(data_lines)), sep="\t", header=None)
    df.columns = [
        "query", "subject", "identity", "alignment_length", "mismatches", "gap_opens",
        "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"
    ]
    df = df.apply(pd.to_numeric, errors="ignore")
    return df, lines

# 筛选函数
def apply_filter_to_all_queries(params):
    all_results = []
    total_hits = 0
    skipped = []
    for qfile in os.listdir(SPLIT_DIR):
        if not qfile.endswith(".txt"):
            continue
        df, _ = load_query_file(os.path.join(SPLIT_DIR, qfile))
        if df.empty:
            skipped.append(qfile)
            continue
        df_filtered = df[
            (df["identity"] > params["identity"]) &
            (df["alignment_length"] >= params["alignment_length"]) &
            (df["mismatches"] <= params["mismatches"]) &
            (df["evalue"] < params["evalue"])
        ]
        df_sorted = df_filtered.sort_values(by=params["sort_column"], ascending=params["ascending"])
        df_topn = df_sorted.groupby("query").head(params["top_n"])
        total_hits += len(df_topn)
        all_results.append(df_topn)
    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
    else:
        combined = pd.DataFrame()
    return combined, total_hits, skipped

# 参数设置（侧栏）
with st.sidebar:
    st.header("📌 筛选参数设置")
    identity_thresh = st.slider("最小 identity (%)", 0.0, 100.0, 90.0)
    alignment_thresh = st.number_input("最小 alignment_length", value=1000)
    mismatches_thresh = st.number_input("最大 mismatches", value=10)
    evalue_thresh = st.number_input("最大 e-value", value=1e-5, format="%.1e")
    top_n = st.number_input("每个 query 保留前 N 个 hits", value=5, min_value=1)
    sort_column = st.selectbox("排序字段", ["identity", "alignment_length", "mismatches", "evalue", "bit_score"])
    ascending = st.radio("排序方式", ["降序", "升序"]) == "升序"

    if st.button("🔄 更新所有 query 的筛选结果"):
        with st.spinner("正在分析所有 query 文件，请稍候..."):
            params = {
                "identity": identity_thresh,
                "alignment_length": alignment_thresh,
                "mismatches": mismatches_thresh,
                "evalue": evalue_thresh,
                "top_n": top_n,
                "sort_column": sort_column,
                "ascending": ascending,
            }
            df_all, total, skipped = apply_filter_to_all_queries(params)
            st.session_state.results_ready = True
            st.session_state.global_result_df = df_all
            st.session_state.total_hits = total
            st.session_state.skipped_queries = skipped
        st.success(f"🎯 已完成筛选：共命中 {total} 条记录")

# 单个 query 展示
selected_file = st.selectbox("选择一个 Query 查看其结果：", query_files)
df_query, raw_lines = load_query_file(os.path.join(SPLIT_DIR, selected_file))

df_q_filtered = df_query[
    (df_query["identity"] > identity_thresh) &
    (df_query["alignment_length"] >= alignment_thresh) &
    (df_query["mismatches"] <= mismatches_thresh) &
    (df_query["evalue"] < evalue_thresh)
]
df_q_sorted = df_q_filtered.sort_values(by=sort_column, ascending=ascending)
df_q_topn = df_q_sorted.groupby("query").head(top_n)

st.subheader(f"🔎 {selected_file} 的筛选结果（Top {top_n}）")
st.dataframe(df_q_topn, use_container_width=True)

# 下载当前 query
tsv_current = df_q_topn.to_csv(index=False, sep="\t")
st.download_button(
    label="📥 下载当前 query 的筛选结果",
    data=tsv_current,
    file_name=f"filtered_{selected_file}",
    mime="text/tab-separated-values"
)

# 所有 query 的合并结果
if st.session_state.results_ready:
    st.subheader("📊 所有 Query 的合并筛选结果")
    st.write(f"✅ 共命中 {st.session_state.total_hits} 条记录")
    st.dataframe(st.session_state.global_result_df.head(50), use_container_width=True)

    all_tsv = st.session_state.global_result_df.to_csv(index=False, sep="\t")
    st.download_button(
        label="📥 下载所有 query 的筛选结果（合并）",
        data=all_tsv,
        file_name="all_filtered_queries.tsv",
        mime="text/tab-separated-values"
    )

    if st.session_state.skipped_queries:
        with st.expander("⚠️ 被跳过的空 query 文件列表"):
            for q in st.session_state.skipped_queries:
                st.text(q)

# 原始注释信息
with st.expander("📝 原始注释信息"):
    for line in raw_lines:
        if line.startswith("#"):
            st.text(line.strip())
