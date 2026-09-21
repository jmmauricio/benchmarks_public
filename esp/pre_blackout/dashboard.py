import sys
from pathlib import Path

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_sld_interactive import build_figure

st.set_page_config(page_title="ESP Pre-Blackout Dashboard", layout="wide")
st.title("ESP Pre-Blackout Dashboard")

tab_sld, tab_eig = st.tabs(["Single-Line Diagram", "Eigenvalues"])

with tab_eig:
    df = pd.read_csv("esp_preblackout_eig.csv")
    col1, col2 = df.columns[1], df.columns[2]

    fig = px.scatter(
        df,
        x=col1,
        y=df[col2] / (2 * np.pi),
        title=f"{col2} vs {col1}",
        hover_data={"Freq": ":.4f", "Damp": ":.4f", "Participation": True},
        labels={col1: col1, col2: col2},
        template="plotly_white",
    )
    fig.update_traces(marker=dict(size=6, opacity=0.7))

    x_min = df[col1].min() * 1.1
    slope_5 = np.sqrt(1 / 0.05**2 - 1)
    slope_10 = np.sqrt(1 / 0.10**2 - 1)
    y_max = df[col2].max() / (2 * np.pi) * 1.05

    for damp_pct, slope, color in [(5, slope_5, "red"), (10, slope_10, "orange")]:
        y_end = slope * abs(x_min) / (2 * np.pi)
        fig.add_shape(
            type="line", x0=0, y0=0, x1=x_min, y1=y_end,
            line=dict(color=color, width=1, dash="dash"),
            layer="below",
        )
        fig.add_annotation(
            x=x_min * 0.95, y=y_end * 0.95,
            text=f"{damp_pct}%", showarrow=False,
            font=dict(color=color, size=11),
        )

    fig.update_yaxes(range=[0, y_max], title_text="Imag (Hz)")
    st.plotly_chart(fig, use_container_width=True)

    with st.expander("Data Table"):
        st.dataframe(df)

with tab_sld:
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        size_by = st.selectbox(
            "Generator size by",
            ["S_n", "P", "Q", "I"],
            format_func=lambda x: {
                "S_n": "Installed capacity (MVA)",
                "P": "Active power P (MW)",
                "Q": "Reactive power |Q| (Mvar)",
                "I": "Current |I| (kA)",
            }[x],
        )
    with c2:
        labels = st.selectbox(
            "Labels",
            ["off", "region", "power"],
            format_func=lambda x: {
                "off": "Off",
                "region": "Region code",
                "power": "Active power (MW)",
            }[x],
        )
    with c3:
        province_fill = st.selectbox(
            "Province fill",
            ["none", "voltage", "active_power", "reactive_power"],
            format_func=lambda x: {
                "none": "None",
                "voltage": "Voltage (pu)",
                "active_power": "Active Power (MW)",
                "reactive_power": "Reactive Power (Mvar)",
            }[x],
        )
    with c4:
        line_color_by = st.selectbox(
            "Line coloring",
            ["none", "active_power", "reactive_power", "apparent_power"],
            format_func=lambda x: {
                "none": "None",
                "active_power": "Active Power (MW)",
                "reactive_power": "Reactive Power (Mvar)",
                "apparent_power": "Apparent Power (MVA)",
            }[x],
        )

    fig_sld = build_figure(
        province_fill=province_fill if province_fill != "none" else None,
        line_color_by=line_color_by if line_color_by != "none" else None,
        size_by=size_by, labels=labels,
    )
    st.plotly_chart(fig_sld, use_container_width=True)
