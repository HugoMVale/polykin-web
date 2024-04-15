import numpy as np
import pandas as pd
import streamlit as st
from polykin.copolymerization.fitting import fit_reactivity_ratios
from streamlit import session_state as state

st.set_page_config(
    page_title="PolyKin | Fitting",
    page_icon="📊",
    layout="wide",
    initial_sidebar_state="expanded")

# Cache stuff
if "init" not in state:
    state.init = True
    state.options_fit_methods = {
        "Nonlinear least squares": "NLLS",
        "Orthogonal distance regression": "ODR"
    }
    state.options_JCR_methods = ['linear', 'exact']

    state.data_default = pd.DataFrame(
        data={
            "f1": [0.186, 0.299, 0.527, 0.600, 0.700, 0.798],
            "F1": [0.196, 0.279, 0.415, 0.473, 0.542, 0.634],
            "scale_f1": [0.01]*6,
            "scale_F1": [0.05]*6
        },
        dtype=np.float64)

st.title('Fit Reactivity Ratios')

with st.sidebar:
    fit_method = st.selectbox("Regression method",
                              state.options_fit_methods.keys(),
                              index=0)

    JCR_method = st.multiselect("Joint confidence region method(s)",
                                state.options_JCR_methods,
                                default=state.options_JCR_methods,
                                format_func=lambda x: x.title())

    cl = st.slider('Confidence level (%)', 80, 95, value=95, step=5)

col = st.columns(2)

with col[0]:
    st.subheader("Experimental Data Input")
    st.text("Enter your copolymer composition data and estimated uncertainty.")

    data = st.data_editor(
        state.data_default,
        column_config={
            "f1": st.column_config.NumberColumn(
                "f₁",
                help="Mole fraction of M1 in the unreacted monomer mixture.",
                min_value=0.,
                max_value=1.),
            "F1": st.column_config.NumberColumn(
                "F₁",
                help="Instantaneous mole fraction of M1 in the copolymer.",
                min_value=0.,
                max_value=1.),
            "scale_f1": st.column_config.NumberColumn(
                "Scale f₁",
                help="Scale factor for f₁. Can be thought of as the absolute standard deviation of f₁.",
                min_value=0.,
                max_value=1.),
            "scale_F1": st.column_config.NumberColumn(
                "Scale F₁",
                help="Scale factor for F₁. Can be thought of as the absolute standard deviation of F₁.",
                min_value=0.,
                max_value=1)
        },
        hide_index=False,
        num_rows='dynamic')

data.dropna(inplace=True)
ndata = data.shape[0]

res = None
if not data.empty:
    try:
        res = fit_reactivity_ratios(
            f1=data['f1'].values,
            F1=data['F1'].values,
            scale_f=data['scale_f1'].values,
            scale_F=data['scale_F1'].values,
            method=state.options_fit_methods[fit_method],
            alpha=(1. - cl/100),
            JCR_method=JCR_method if ndata > 2 else []
        )
    except ValueError as e:
        col[0].error(
            "Oops, we have a problem. Check the data and/or turn off the exact JCR.")
        with col[0].expander("Error message"):
            st.write(e.args)
else:
    col[0].warning("Try adding some data.")

if res:
    with col[0]:

        st.subheader("Results")

        table_results = {
            "Parameter": ["r1", "r2"],
            "Point estimate": [res.r1, res.r2],
            "Standard error": [res.se_r1, res.se_r2],
            f"Confidence interval {cl}%": [res.ci_r1, res.ci_r2]
        }
        st.dataframe(table_results, hide_index=True)

        if ndata > 2:
            st.write("Covariance matrix")
            st.dataframe(pd.DataFrame(res.cov, columns=["r1", "r2"],
                                      index=["r1", "r2"]))

    with col[1]:
        Mayo = res.Mayo
        if Mayo and Mayo[0]:
            st.pyplot(Mayo[0])

        JCR = res.JCR
        if JCR and JCR[0]:
            st.pyplot(JCR[0])
