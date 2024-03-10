import os
import scanpy as sc
import scTenifoldXct as sct
import streamlit as st
from tempfile import NamedTemporaryFile
from datetime import datetime
from utils import *


KEYS = ["data", "data_obs", "sender", "receiver", "rebuild_on"]
DATE = datetime.now().strftime("%m%d%y")


# load data
st.header("🤖 scTenifold for single-sample analysis")
st.info("This is a beta version UI for [scTenifoldXct](https://github.com/cailab-tamu/scTenifoldXct) project")

st.write("Upload your h5ad file")
uploaded_file = st.file_uploader("File upload", type="h5ad")
if uploaded_file is not None:
    with NamedTemporaryFile(dir=".", suffix=".h5ad") as f:
        f.write(uploaded_file.getbuffer())
        # st.write(f.name) # file path of temp.h5ad
        data = sc.read_h5ad(f.name)
        if "data" not in st.session_state:
            st.session_state["data"] = data
        if "data_obs" not in st.session_state:
            st.session_state["data_obs"] = data.obs.columns.to_list()

    # cell type column
    option_obs = st.selectbox(
    "Select column name containing cell types",
    st.session_state["data_obs"],
    index=None,
    placeholder="Select column...",
    )
    if option_obs is not None :
        st.session_state["option_obs"] = option_obs
        # st.write('You selected:', option_obs)
        # sender and receiver cells
        cell_types = st.session_state["data"].obs[st.session_state["option_obs"]].unique().to_list()
        if "cell_types" not in st.session_state:
            st.session_state["cell_types"] = cell_types

        sender = st.selectbox(
        "Select sender cell type",
        st.session_state["cell_types"],
        index=None,
        placeholder="Select sender cell type...",
        )
        if sender is not None:
            st.session_state["sender"] = sender
        # st.write('You selected:', sender)

        receiver = st.selectbox(
        "Select receiver cell type",
        st.session_state["cell_types"],
        index=None,
        placeholder="Select receiver cell type...",
        )
        if receiver is not None:
            st.session_state["receiver"] = receiver
        # st.write('You selected:', receiver)

        if sender is not None and receiver is not None:
            # build GRNs
            rebuild_on = st.selectbox(
            "Select whether to build GRNs",
            [True, False],
            index=None,
            placeholder="Select False for rerun...",
            )
            if rebuild_on is not None:
                st.session_state["rebuild_on"] = rebuild_on
                n_cpu = st.slider(
                "Choose CPUs to use",
                min_value=1,
                max_value=os.cpu_count(),
                step=1,
                )
                if n_cpu not in st.session_state:
                    st.session_state["n_cpu"]=n_cpu

        if check_values_not_none(st.session_state, KEYS):
            # run
            run_ = st.button("Click and run scTenifoldXct", type="primary")
            complete_ = False
            if run_:
                with st.spinner("Running..."):
                    xct = sct.scTenifoldXct(data = st.session_state["data"], # an AnnData 
                                            source_celltype = st.session_state["sender"], # sender cell type
                                            target_celltype = st.session_state["receiver"], # receiver cell type
                                            obs_label = st.session_state["option_obs"], # colname in adata.obs indicating cell types
                                            rebuild_GRN = st.session_state["rebuild_on"], # whether to build GRNs
                                            GRN_file_dir = "GRNs",  # folder path to GRNs
                                            verbose = True, # whether to verbose the processing
                                            n_cpus = st.session_state["n_cpu"]) # CPU multiprocessing, -1 to use all
                    result_df = get_xct_result(xct)
                    complete_ = True
                st.success("Done!")
            if complete_:  
                result_csv = convert_df(result_df)
                st.download_button(
                    label="Download result as CSV",
                    data=result_csv,
                    file_name=f"xct_result_{DATE}.csv",
                    mime="text/csv",
                )

with st.expander("Implementing Tips"):
    st.write('''
    - Training may take several minutes
    - Check our [paper](https://doi.org/10.1016/j.cels.2023.01.004) and [repository](https://github.com/cailab-tamu/scTenifoldXct) for more information
    ''')

st.divider()
st.write('''
(c) 2024 CaiLab, Texas A&M University
''')
# st.write("session_state:", st.session_state) # debug

# streamlit run xctapp.py
                