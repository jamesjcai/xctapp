import os
import scanpy as sc
import scTenifoldXct as sct
import streamlit as st
from tempfile import NamedTemporaryFile
from datetime import datetime
from utils import *


DATE = datetime.now().strftime("%m%d%y")

st.header("🤖 scTenifold for two-sample analysis")
st.info("This is a beta version UI for [scTenifoldXct](https://github.com/cailab-tamu/scTenifoldXct) project")

def add_item(col, key):
    col.write(f"Upload h5ad file for {key}")
    uploaded_file = col.file_uploader("File upload", type="h5ad", key=key)
    if uploaded_file is not None:
        with NamedTemporaryFile(dir=".", suffix=".h5ad") as f:
            f.write(uploaded_file.getbuffer())
            data = sc.read_h5ad(f.name)
            if f"{key}_data" not in st.session_state:
                st.session_state[f"{key}_data"] = data
            if f"{key}_data_obs" not in st.session_state:
                st.session_state[f"{key}_data_obs"] = data.obs.columns.to_list()

        # cell type column
        option_obs = st.selectbox(
        "Select column name containing cell types",
        st.session_state[f"{key}_data_obs"],
        index=None,
        placeholder="Select column...",
        key=f"{key}_opt1",
        )
        if option_obs is not None :
            # sender and receiver cells
            cell_types = st.session_state[f"{key}_data"].obs[st.session_state[f"{key}_opt1"]].unique().to_list()
            if f"{key}_cell_types" not in st.session_state:
                st.session_state[f"{key}_cell_types"] = cell_types

            sender = st.selectbox(
            "Select sender cell type",
            st.session_state[f"{key}_cell_types"],
            index=None,
            placeholder="Select sender cell type...",
            key=f"{key}_opt2",
            )

            receiver = st.selectbox(
            "Select receiver cell type",
            st.session_state[f"{key}_cell_types"],
            index=None,
            placeholder="Select receiver cell type...",
            key=f"{key}_opt3",
            )

            if sender is not None and receiver is not None:
                # build GRNs
                rebuild_on = st.selectbox(
                "Select whether to build GRNs",
                [True, False],
                index=None,
                placeholder="Select False for rerun...",
                key=f"{key}_opt4",
                )
                if rebuild_on is not None:
                    n_cpu = st.slider(
                    "Choose CPUs to use",
                    min_value=1,
                    max_value=os.cpu_count(),
                    step=1,
                    key=f"{key}_opt5",
                    )
                    if f"{key}_opt5" in st.session_state:
                        create_ = st.button(f"Create {key}", key=f"{key}_create")
                        complete_ = False
                        if create_:
                            xct = sct.scTenifoldXct(data = st.session_state[f"{key}_data"], # an AnnData 
                                                    source_celltype = st.session_state[f"{key}_opt2"], # sender cell type
                                                    target_celltype = st.session_state[f"{key}_opt3"], # receiver cell type
                                                    obs_label = st.session_state[f"{key}_opt1"], # colname in adata.obs indicating cell types
                                                    rebuild_GRN = st.session_state[f"{key}_opt4"], # whether to build GRNs
                                                    GRN_file_dir = f"GRNs_{key}",  # folder path to GRNs # TODO
                                                    verbose = True, # whether to verbose the processing
                                                    n_cpus = st.session_state[f"{key}_opt5"]) # CPU multiprocessing, -1 to use all
                            st.session_state[f"{key}_xct_object"] = xct
                            complete_ = True
                            if f"{key}_complete" not in st.session_state:
                                st.session_state[f"{key}_complete"] = complete_
                            st.success(f"{key} created!", icon="✅")

cols = st.columns(2)
for i, col in enumerate(cols):
    with col:
        add_item(col=col, key=f"item{i+1}")

if "item1_complete" in st.session_state and "item2_complete" in st.session_state:
    complete_ = False
    run_ = st.button("Click and run scTenifoldXct diff", key="run", type="primary")
    if run_:
        with st.spinner("Running..."):
            xcts_pairs_diff = get_merge_xct_result(st.session_state[f"item1_xct_object"], st.session_state[f"item2_xct_object"])
            complete_ = True
        st.success(f"Done!")
    if complete_:  
        xcts_pairs_diff = convert_df(xcts_pairs_diff)
        st.download_button(
            label="Download result as CSV",
            data=xcts_pairs_diff,
            file_name=f"xct_diff_result_{DATE}.csv",
            mime="text/csv",
        )

with st.expander("Implementing Tips"):
    st.write('''
    - Training for two-sample takes longer time than simple-sample
    - Check our [paper](https://doi.org/10.1016/j.cels.2023.01.004) and [repository](https://github.com/cailab-tamu/scTenifoldXct) for more information
    ''')

st.divider()
st.write('''
(c) 2024 CaiLab, Texas A&M University
''')
# st.write("session_state:", st.session_state) # debug
