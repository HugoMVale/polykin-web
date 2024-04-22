from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd
import streamlit as st

# %% Parameters

options_fit_methods = {
    "Nonlinear least squares": "NLLS",
    "Orthogonal distance regression": "ODR"
}

options_JCR_methods = ['linear', 'exact']

latex_symbols = {
    "Ff": "$F_1(f_1)$",
    "fx": "$f_1(x)$",
    "Fx": "$F_1(x)$",
    "r1": "$r_1$",
    "r2": "$r_2$"
}

# %% Data stuff

default_tables = {
    "Ff": pd.DataFrame(data={"f1": [0.186, 0.299, 0.527, 0.600, 0.700, 0.798],
                             "F1": [0.196, 0.279, 0.415, 0.473, 0.542, 0.634],
                             "scale_f1": [0.01]*6,
                             "scale_F1": [0.05]*6},
                       dtype=np.float64),
    "fx": pd.DataFrame(data={"x": [0.], "f1": [0.5], "scale_f1": [1.]},
                       dtype=np.float64),
    "Fx": pd.DataFrame(data={"x": [0.], "F1": [0.5], "scale_F1": [1.]},
                       dtype=np.float64),

}

column_config = {
    "Ff": {
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
            max_value=1.)
    },
    "fx": {
        "x": st.column_config.NumberColumn(
            "x",
            help="Total monomer molar conversion.",
            min_value=0.,
            max_value=1.),
        "f1": st.column_config.NumberColumn(
            "f₁",
            help="Mole fraction of M1 in the unreacted monomer mixture.",
            min_value=0.,
            max_value=1.),
        "scale_f1": st.column_config.NumberColumn(
            "Scale f₁",
            help="Scale factor for f₁. Can be thought of as the absolute standard deviation of f₁.",
            min_value=0.,
            max_value=1.)
    },
    "Fx": {
        "x": st.column_config.NumberColumn(
            "x",
            help="Total monomer molar conversion.",
            min_value=0.,
            max_value=1.),
        "F1": st.column_config.NumberColumn(
            "F₁",
            help="Cumulative mole fraction of M1 in the copolymer.",
            min_value=0.,
            max_value=1.),
        "scale_F1": st.column_config.NumberColumn(
            "Scale F₁",
            help="Scale factor for F₁. Can be thought of as the absolute standard deviation of F₁.",
            min_value=0.,
            max_value=1.)
    }
}


@dataclass
class DataSet():
    name: str = ""
    weight: float = 1.
    f10: float = 1.
    data: Optional[pd.DataFrame] = None


def update_data(datasets: list[DataSet],
                n: int
                ) -> None:
    ndata = len(datasets)
    if n < ndata:
        del datasets[n:]
    elif n > ndata:
        nadd = n - ndata
        datasets.extend([DataSet()]*nadd)
    else:
        pass
