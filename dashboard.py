import sqlite3
from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

import analysis

DB_PATH = Path(__file__).resolve().parent / "cell-count.db"

RESPONSE_COLORS = {"yes": "#2a78d6", "no": "#eb6834"}
CATEGORY_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]

st.set_page_config(page_title="Loblaw Bio", layout="wide")


@st.cache_resource
def get_connection() -> sqlite3.Connection:
    if not DB_PATH.exists():
        st.error(f"Database not found at {DB_PATH}. Run `python load_data.py` first.")
        st.stop()
    conn = sqlite3.connect(DB_PATH, check_same_thread=False)
    conn.row_factory = sqlite3.Row
    return conn


@st.cache_data
def load_relative_frequencies() -> pd.DataFrame:
    return pd.DataFrame(analysis.relative_frequencies(get_connection()))


@st.cache_data
def load_responder_frequencies() -> pd.DataFrame:
    return pd.DataFrame(analysis.melanoma_miraclib_pbmc_frequencies(get_connection()))


@st.cache_data
def load_responder_significance() -> pd.DataFrame:
    return pd.DataFrame(analysis.responder_significance(get_connection()))


@st.cache_data
def load_baseline_samples() -> pd.DataFrame:
    return pd.DataFrame(analysis.baseline_melanoma_miraclib_pbmc_samples(get_connection()))


@st.cache_data
def load_baseline_summary() -> pd.DataFrame:
    return pd.DataFrame(analysis.baseline_subset_summary(get_connection()))


st.title("Immune Cell Population Dashboard")

tab1, tab2, tab3 = st.tabs(
    [
        "Initial analysis",
        "Statistical analysis",
        "Data subset analysis",
    ]
)

with tab1:
    st.subheader("Relative frequency of each cell population per sample")
    df2 = load_relative_frequencies()

    col1, col2 = st.columns([1, 2])
    with col1:
        populations = st.multiselect(
            "Population", sorted(df2["population"].unique()), key="p2_pop"
        )
    with col2:
        sample_search = st.text_input("Search sample ID", key="p2_sample")

    filtered = df2
    if populations:
        filtered = filtered[filtered["population"].isin(populations)]
    if sample_search:
        filtered = filtered[filtered["sample"].str.contains(sample_search, case=False)]

    st.dataframe(filtered, use_container_width=True, height=400)
    st.download_button(
        "Download full table as CSV",
        df2.to_csv(index=False),
        file_name="part_2_relative_frequencies.csv",
        mime="text/csv",
    )
    st.caption(f"Showing {len(filtered):,} of {len(df2):,} rows")

with tab2:
    st.subheader("Melanoma / miraclib / PBMC samples: relative frequency by responder status")

    df3a = load_responder_frequencies()
    fig = px.box(
        df3a,
        x="response",
        y="percentage",
        color="response",
        facet_col="population",
        category_orders={"response": ["yes", "no"]},
        color_discrete_map=RESPONSE_COLORS,
        labels={"percentage": "Relative frequency (%)", "response": "Responder"},
    )
    fig.update_xaxes(matches=None)
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    st.plotly_chart(fig, use_container_width=True)
    st.download_button(
            "Download data as CSV",
            df3a.to_csv(index=False),
            file_name="part_3_melanoma_miraclib_pbmc_frequencies_frequencies.csv",
            mime="text/csv",
        )

    st.subheader("Statistical significance (Mann-Whitney U, BH-adjusted)")
    df3b = load_responder_significance()
    st.dataframe(
        df3b.style.apply(
            lambda row: ["background-color: rgba(42,120,214,0.15)"] * len(row)
            if row["significant"]
            else [""] * len(row),
            axis=1,
        ),
        use_container_width=True,
    )
    st.caption(
        "No population remains significant after correcting for testing 5 populations "
        "at once (Benjamini-Hochberg FDR, alpha = 0.05)."
        if not df3b["significant"].any()
        else "Highlighted rows are significant after FDR correction (alpha = 0.05)."
    )

with tab3:
    st.subheader("Baseline (day 0) melanoma / miraclib / PBMC sample list")

    df4a = load_baseline_samples()
    
    col1, col2, col3 = st.columns(3)
    with col1:
        project_filter = st.multiselect("Project", sorted(df4a["project"].unique()))
    with col2:
        response_filter = st.multiselect("Response", sorted(df4a["response"].unique()))
    with col3:
        sex_filter = st.multiselect("Sex", sorted(df4a["sex"].unique()))

    filtered4 = df4a
    if project_filter:
        filtered4 = filtered4[filtered4["project"].isin(project_filter)]
    if response_filter:
        filtered4 = filtered4[filtered4["response"].isin(response_filter)]
    if sex_filter:
        filtered4 = filtered4[filtered4["sex"].isin(sex_filter)]

    st.dataframe(filtered4, use_container_width=True, height=400)
    st.download_button(
        "Download full table as CSV",
        df4a.to_csv(index=False),
        file_name="part_4a_baseline_samples.csv",
        mime="text/csv",
    )
    st.caption(f"Showing {len(filtered4):,} of {len(df4a):,} rows")

    st.subheader("Baseline sample subsets")

    df4b = load_baseline_summary()
    col1, col2, col3 = st.columns(3)
    for col, category, title in zip(
        (col1, col2, col3),
        ("project", "response", "sex"),
        ("Samples per project", "Subjects per response", "Subjects per sex"),
    ):
        with col:
            subset = df4b[df4b["category"] == category]
            fig = px.bar(
                subset,
                x="value",
                y="count",
                color="value",
                color_discrete_sequence=CATEGORY_COLORS,
                text="count",
            )
            fig.update_layout(showlegend=False, title=title, xaxis_title=None)
            st.plotly_chart(fig, use_container_width=True)
