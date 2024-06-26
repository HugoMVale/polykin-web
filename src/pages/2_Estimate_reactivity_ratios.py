import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
from matplotlib.figure import Figure
from polykin.copolymerization import (CopoDataset_Ff, CopoDataset_fx,
                                      CopoDataset_Fx, fit_copo_data)
from polykin.copolymerization.binary import (inst_copolymer_binary,
                                             monomer_drift_binary)
from streamlit import session_state as state

import utils

st.set_page_config(
    page_title="Copo Fitting | polykin",
    page_icon="📊",
    layout="wide",
    initial_sidebar_state="expanded")

if "page" not in state or state.page != 2:
    state.page = 2
    state.data = {"Ff": [], "fx": [], "Fx": []}
    state.counter = 0

# %% Data and Settings - left column

column = st.columns((0.8, 1, 1), gap='large')

with column[0]:

    st.header("Experimental Data")

    with st.expander("Enter the type and number of datasets"):
        number_sets = {}
        for kind in utils.default_tables.keys():
            symbol = utils.latex_symbols[kind]
            number_sets[kind] = st.number_input(
                f"Number of {symbol} datasets",
                min_value=0, max_value=5, value=0, step=1,
                key=kind + "_number")

with st.sidebar:

    options_fit_methods = [*utils.options_fit_methods.keys()]
    if (number_sets['fx'] + number_sets['Fx']) > 0:
        options_fit_methods.pop(-1)

    fit_method = st.selectbox("Regression method",
                              options_fit_methods,
                              index=0)
    ratios_guess = [1., 1.]
    for idx, name in enumerate(["r1", "r2"]):
        key = name + "_guess"
        symbol = utils.latex_symbols[name]
        ratios_guess[idx] = st.number_input(
            f"Initial guess {symbol}",
            min_value=0., max_value=100., value=1., step=0.1,
            key=name+"_guess")

    confidence_level = st.slider('Confidence level (%)', min_value=80,
                                 max_value=95, value=95, step=5, key="CL")
    JCR_linear = st.toggle("Compute linear JCR", value=True)
    JCR_exact = st.toggle("Compute exact JCR", value=False)


# Ajust dataset number
for kind in number_sets.keys():
    datasets = state.data[kind]
    utils.update_data(datasets, int(number_sets[kind]))


with column[0]:

    for kind in number_sets.keys():

        if number_sets[kind]:

            st.subheader(f"{utils.latex_symbols[kind]} datasets")
            datasets = state.data[kind]

            for idx, ds in enumerate(datasets):

                basekey = f"{kind}-{idx+1}"

                ds.name = st.text_input(
                    "Dataset name", value=basekey,
                    key=basekey+"_name")

                ds.weight = st.slider(
                    "Relative dataset weight (%)",
                    min_value=0, max_value=100, value=100, step=10,
                    key=basekey+"_weight")

                if kind != "Ff":
                    ds.f10 = st.number_input(
                        r"$f_{1}(0)$",
                        min_value=0., max_value=1., value=0.5,
                        help="Initial comonomer composition",
                        key=basekey+"_f10")

                ds.data = st.data_editor(
                    utils.default_tables[kind],
                    column_config=utils.column_config[kind],
                    hide_index=True,
                    num_rows='dynamic',
                    key=basekey+"_data")
                # st.write(ds)
                st.divider()


# %% Plots - middle column

@st.cache_data
def monomer_drift_binary_cached(*args):
    state.counter += 1
    return monomer_drift_binary(*args)


@ st.cache_data
def make_plot_Ff(datasets: list[utils.DataSet],
                 r: tuple[float, float]) -> Figure:
    """Generate F(f) plot."""
    fig, ax = plt.subplots()
    ax.set_xlabel(r"$f_1$")
    ax.set_ylabel(r"$F_1$")
    ax.set_xlim(0., 1.)
    ax.set_ylim(0., 1.)
    for dataset in datasets:
        if dataset.data is not None and not dataset.data.empty:
            f1 = dataset.data['f1']
            F1 = dataset.data['F1']
            ax.scatter(f1, F1, label=dataset.name)
    f1 = np.linspace(0., 1., 200)
    F1 = inst_copolymer_binary(f1, *r)
    ax.plot(f1, F1)
    ax.legend(loc='best')
    return fig


@ st.cache_data
def make_plot_fx(datasets: list[utils.DataSet],
                 r: tuple[float, float]
                 ) -> Figure:
    """Generate f(x) plot."""
    fig, ax = plt.subplots()
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$f_1$")
    ax.set_xlim(0., 1.)
    ax.set_ylim(0., 1.)
    xs = np.linspace(0., 1., 500)
    for dataset in datasets:
        if dataset.data is not None and not dataset.data.empty:
            # exp data
            x = dataset.data['x']
            f1 = dataset.data['f1']
            ax.scatter(x, f1, label=dataset.name)
            # simulated
            f1 = monomer_drift_binary_cached(dataset.f10, xs, *r)
            ax.plot(xs, f1)
    ax.legend(loc='best')
    return fig


@ st.cache_data
def make_plot_Fx(datasets: list[utils.DataSet],
                 r: tuple[float, float]
                 ) -> Figure:
    """Generate F(x) plot."""
    fig, ax = plt.subplots()
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$F_1$")
    ax.set_xlim(0., 1.)
    ax.set_ylim(0., 1.)
    xs = np.linspace(0., 1., 500)
    for dataset in datasets:
        if dataset.data is not None and not dataset.data.empty:
            # exp data
            x = dataset.data['x']
            F1 = dataset.data['F1']
            ax.scatter(x, F1, label=dataset.name)
            # simulated
            f10 = dataset.f10
            f1 = monomer_drift_binary_cached(f10, xs, *r)
            F1 = f1 + (f10 - f1)/(xs + 1e-10)
            ax.plot(xs, F1)
    ax.legend(loc='best')
    return fig


make_plot = {"Ff": make_plot_Ff, "fx": make_plot_fx, "Fx": make_plot_Fx}


def update_plots(r: tuple[float, float]):
    """Update plots."""
    for kind in ["Ff", "fx", "Fx"]:
        if number_sets[kind] > 0:
            st.pyplot(make_plot[kind](state.data[kind], r))


if sum(number_sets.values()) > 0:
    with column[1]:
        st.header("Data Overview")
        plot_zone = st.empty()
        with plot_zone.container():
            update_plots(tuple(ratios_guess))  # type: ignore


# %% Fit


@st.cache_data
def check_data(data: dict[str, list[utils.DataSet]]) -> tuple[list, int]:
    """Check which datasets have missing values and return total number
    of valid rows."""
    datasets_with_null = []
    ndata = 0
    for datasets in data.values():
        for dataset in datasets:
            if dataset.data is not None:
                nulls = dataset.data.isnull()
                found_null = np.count_nonzero(nulls) > 0
                if found_null:
                    datasets_with_null.append(dataset.name)
                ndata += dataset.data.shape[0] - nulls.any(axis=1).sum()
    return (datasets_with_null, ndata)


run_fit = False
if sum(number_sets.values()) > 0:
    with column[2]:
        st.header("Results")
        info_zone = st.empty()

        datasets_with_null, ndata = check_data(state.data)
        if datasets_with_null:
            with info_zone.container():
                st.warning(
                    "The following datasets have missing data: " +
                    ", ".join(datasets_with_null))
        sufficient_data = ndata >= 2
        run_fit = st.button("Run Fit ▷", type="primary",
                            disabled=not (sufficient_data))

if run_fit:

    # Pack data
    data_Ff = []
    data_fx = []
    data_Fx = []
    if state.data['Ff']:
        for ds in state.data['Ff']:
            ds.data.dropna(inplace=True)
            dsi = CopoDataset_Ff(
                name=ds.name,
                f1=ds.data['f1'].values,
                F1=ds.data['F1'].values,
                scale_f1=ds.data['scale_f1'].values,
                scale_F1=ds.data['scale_F1'].values,
                weight=ds.weight)
            data_Ff.append(dsi)

    if state.data['fx']:
        for ds in state.data['fx']:
            ds.data.dropna(inplace=True)
            dsi = CopoDataset_fx(
                name=ds.name,
                f10=ds.f10,
                x=ds.data['x'].values,
                f1=ds.data['f1'].values,
                scale_f1=ds.data['scale_f1'].values,
                weight=ds.weight)
            data_fx.append(dsi)

    if state.data['Fx']:
        for ds in state.data['Fx']:
            ds.data.dropna(inplace=True)
            dsi = CopoDataset_Fx(
                name=ds.name,
                f10=ds.f10,
                x=ds.data['x'].values,
                F1=ds.data['F1'].values,
                scale_F1=ds.data['scale_F1'].values,
                weight=ds.weight)
            data_Fx.append(dsi)

    res = None
    try:
        res = fit_copo_data(data_Ff=data_Ff,
                            data_fx=data_fx,
                            data_Fx=data_Fx,
                            r_guess=tuple(ratios_guess),
                            alpha=(1. - confidence_level/100),
                            plot_data=False,
                            JCR_linear=JCR_linear,
                            JCR_exact=JCR_exact
                            )
    except ValueError as e:
        with info_zone.container():
            st.error(
                "Oopsie, we have a problem... Check the data and/or turn off the exact JCR.")
            with st.expander("Error message"):
                st.write(e.args)

    if res:

        table_results = {
            "Parameter": ["r₁", "r₂"],
            "Point estimate": [res.r1, res.r2],
            "Standard error": [res.se_r1, res.se_r2],
            f"Confidence interval {confidence_level}%": [res.ci_r1, res.ci_r2]
        }

        with plot_zone.container():
            update_plots((res.r1, res.r2))

        with column[2]:

            st.write("Parameter estimates and confidence intervals")
            st.dataframe(table_results, hide_index=True)

            if ndata > 2:
                st.write("Covariance matrix")
                st.dataframe(pd.DataFrame(res.cov, columns=["r₁", "r₂"],
                                          index=["r₁", "r₂"]))

            JCR = res.plots['JCR'][0]
            if JCR is not None:
                st.write("Joint confidence region")
                st.pyplot(JCR)

# %% Save state
# state.ratios_guess_last = ratios_guess
