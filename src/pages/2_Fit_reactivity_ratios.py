import numpy as np
import pandas as pd
import streamlit as st
from polykin.copolymerization.fitting import fit_reactivity_ratios

#
options_fit_methods = {
    "Nonlinear least squares": "NLLS",
    "Orthogonal distance regression": "ODR"
}

options_JCR_methods = ['linear', 'exact']

st.title('Fit Reactivity Ratios')

st.subheader("Experimental Data Input")
st.text("Enter your copolymer composition data and estimated uncertainty.")

fit_method = st.sidebar.selectbox("Regression method",
                                  options_fit_methods.keys(),
                                  index=0)

JCR_method = st.sidebar.multiselect("Joint confidence region method(s)",
                                    options_JCR_methods,
                                    default=options_JCR_methods,
                                    format_func=lambda x: x.title()
                                    )
cl = st.sidebar.slider('Confidence level (%)', 80, 95, value=95, step=5)


# Default composition data
f1_default = [0.186, 0.299, 0.527, 0.600, 0.700, 0.798]
F1_default = [0.196, 0.279, 0.415, 0.473, 0.542, 0.634]

data_default = pd.DataFrame(
    data={
        "f1": f1_default,
        "F1": F1_default,
        "scale_f1": [0.01]*len(f1_default),
        "scale_F1": [0.05]*len(F1_default)
    },
    dtype=np.float64)

data = st.data_editor(data_default,
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

st.subheader("Results")

if not data.empty:
    res = fit_reactivity_ratios(
        f1=data['f1'].values,
        F1=data['F1'].values,
        scale_f=data['scale_f1'].values,
        scale_F=data['scale_F1'].values,
        method=options_fit_methods[fit_method],
        alpha=(1. - cl/100),
        JCR_method=JCR_method if ndata > 2 else []
    )

else:
    res = None
    st.write("Try adding some data.")

    # Results
if res:
    table_results = {
        "Parameter": ["r1", "r2"],
        "Point estimate": [res.r1, res.r2],
        "Standard error": [res.se_r1, res.se_r2],
        f"Confidence interval {cl}%": [res.ci_r1, res.ci_r2]
    }
    st.dataframe(table_results, hide_index=True)

    st.write("Covariance matrix")
    st.dataframe(pd.DataFrame(res.cov, columns=[
                 "r1", "r2"], index=["r1", "r2"]))

    Mayo = res.Mayo
    JCR = res.JCR

    if Mayo and Mayo[0]:
        st.pyplot(Mayo[0])

    if JCR and JCR[0]:
        st.pyplot(JCR[0])
