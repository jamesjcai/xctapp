import streamlit as st
import scTenifoldXct as sct


def check_values_not_none(dictionary, keys):
    return all(key in dictionary and dictionary[key] is not None for key in keys)


@st.cache_data
def convert_df(df):
    return df.to_csv().encode('utf-8')


def get_xct_result(xct):
    emb = xct.get_embeds(train = True) 
    result_df = xct.null_test()
    return result_df


def get_merge_xct_result(xct1, xct2):
    xct_merge = sct.merge_scTenifoldXct(xct1, xct2)
    emb = xct_merge.get_embeds(train = True)
    xct_merge.nn_aligned_diff(emb) 
    xcts_pairs_diff = xct_merge.chi2_diff_test()
    return xcts_pairs_diff
