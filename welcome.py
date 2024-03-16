import streamlit as st

st.set_page_config(
    page_title="About",
    page_icon="👋",
)

st.write("# Welcome to scTenifold land! 👋")

st.sidebar.success("Select a mode above.")

st.markdown(
    """
    ScTenifoldXct is an open-source framework designed specifically for inferring cell-cell interactions.\\
    **👈 Select a mode from the sidebar** 
    ### Want to learn more?
    - Check out our [paper](https://doi.org/10.1016/j.cels.2023.01.004)
    - Jump into our [repository](https://github.com/cailab-tamu/scTenifoldXct)
"""
)

st.divider()
st.write('''
(c) 2024 CaiLab, Texas A&M University
''')
