from __future__ import annotations

import sqlite3
from pathlib import Path

import pandas as pd
import streamlit as st
import plotly.express as px

from analysis import (
    compute_relative_frequencies,
    compare_responders,
    baseline_summary,
)


@st.cache_resource
def load_database(db_path: str | Path) -> Dict[str, pd.DataFrame] | Dict[str, Any]:
    """Load and prepare analysis data from the SQLite database.

    This function reads the provided SQLite database and returns a
    dictionary of data structures used throughout the dashboard.  It
    performs the following computations:

    - **Relative frequency summary**: uses :func:`compute_relative_frequencies`
      to calculate the total cell count per sample, the count per immune
      population and the relative frequency (percentage).  This produces
      a long-format DataFrame with one row per (sample, population)
      combination.
    - **Responder comparison**: calls :func:`compare_responders` to
      perform statistical comparisons of relative frequencies between
      responders and non‑responders among melanoma patients receiving
      miraclib.
    - **Baseline statistics**: obtains baseline summaries via
      :func:`baseline_summary`.
    - **Raw data table**: constructs a wide-format DataFrame with one
      row per sample.  The table includes sample metadata (project,
      subject, condition, age, sex, treatment, response, sample type,
      time from treatment start) and one column for each immune cell
      population count.  This structure makes it easy to filter the
      underlying data by metadata fields in the dashboard.

    Parameters
    ----------
    db_path : str or Path
        Path to the SQLite database file.

    Returns
    -------
    dict
        A dictionary with four keys:

        ``"rel_freq"`` : pandas.DataFrame
            Relative frequency summary table.

        ``"stats_df"`` : pandas.DataFrame
            Statistical comparison of responders vs non‑responders.

        ``"baseline"`` : dict
            Baseline summary statistics.

        ``"raw_data"`` : pandas.DataFrame
            Wide-format table containing sample metadata and cell counts.
    """
    # Compute relative frequencies
    rel_freq = compute_relative_frequencies(db_path)
    # Statistical comparison between responders vs non responders
    stats_df = compare_responders(db_path)
    # Baseline summary
    baseline = baseline_summary(db_path)

    # Construct raw_data: join samples with cell_counts and pivot
    conn = sqlite3.connect(db_path)
    # Pull a long-format dataset: one row per measurement per sample
    df_raw = pd.read_sql_query(
        """
        SELECT
            s.id AS sample_id,
            s.sample,
            s.project,
            s.subject,
            s.condition,
            s.age,
            s.sex,
            s.treatment,
            s.response,
            s.sample_type,
            s.time_from_treatment_start,
            c.population,
            c.count
        FROM samples AS s
        LEFT JOIN cell_counts AS c
            ON s.id = c.sample_id
        """,
        conn,
    )
    conn.close()

    # Pivot the long-format counts into a wide table: one row per sample
    if not df_raw.empty:
        pivot = df_raw.pivot_table(
            index=[
                "sample_id",
                "sample",
                "project",
                "subject",
                "condition",
                "age",
                "sex",
                "treatment",
                "response",
                "sample_type",
                "time_from_treatment_start",
            ],
            columns="population",
            values="count",
            fill_value=0,
        ).reset_index()
        # Flatten multi-index columns if present
        pivot.columns = [
            col if not isinstance(col, tuple) else col[1] for col in pivot.columns
        ]
        raw_data = pivot
    else:
        # No data available; create an empty DataFrame with the same
        # metadata columns as df_raw (without population and count).
        metadata_cols = [
            "sample_id",
            "sample",
            "project",
            "subject",
            "condition",
            "age",
            "sex",
            "treatment",
            "response",
            "sample_type",
            "time_from_treatment_start",
        ]
        raw_data = df_raw[metadata_cols].drop_duplicates().copy()

    return {
        "rel_freq": rel_freq,
        "stats_df": stats_df,
        "baseline": baseline,
        "raw_data": raw_data,
    }


def main():
    st.set_page_config(
        page_title="Immune Cell Population Dashboard",
        layout="wide",
    )
    st.title("Immune Cell Population Analysis")
    st.markdown(
        """
        Use the sidebar to navigate between different analyses.  This dashboard
        summarises cell counts for various immune populations across
        samples, compares relative frequencies between responders and
        non‑responders in a melanoma trial, and explores baseline
        characteristics of the dataset.
        """
    )
    # Sidebar
    st.sidebar.header("Navigation")
    page = st.sidebar.radio(
        "Select a view:",
        options=["Data Overview", "Responder Analysis", "Baseline Summary"],
    )
    # Database path selection
    st.sidebar.subheader("Database")
    default_db = Path(__file__).resolve().parent / "immune_data.db"
    db_path = st.sidebar.text_input("Path to SQLite database", value=str(default_db))
    if not Path(db_path).exists():
        st.sidebar.error(f"Database not found at {db_path}. Please generate it first.")
        return
    # Load data
    data = load_database(db_path)
    rel_freq = data["rel_freq"]
    stats_df = data["stats_df"]
    baseline = data["baseline"]
    # Data Overview
    if page == "Data Overview":
        st.header("Relative Frequencies per Sample")
        st.write(
            "This table shows, for each sample, the total cell count, population count and relative frequency (percentage) of each immune cell population."
        )
        # Option to filter by sample or population
        samples = sorted(rel_freq["sample"].unique().tolist())
        selected_samples = st.multiselect(
            "Filter by sample (optional)", options=samples, default=[]
        )
        populations = sorted(rel_freq["population"].unique().tolist())
        selected_populations = st.multiselect(
            "Filter by population (optional)", options=populations, default=populations
        )
        df_display = rel_freq.copy()
        if selected_samples:
            df_display = df_display[df_display["sample"].isin(selected_samples)]
        if selected_populations:
            df_display = df_display[df_display["population"].isin(selected_populations)]
        st.dataframe(
            df_display.reset_index(drop=True),
            use_container_width=True,
            height=500,
        )

        # Raw data view with filtering
        st.subheader("Raw Data (samples and counts)")
        st.write(
            "View the full sample-level data with cell counts and apply filters across metadata fields."
        )
        # Load raw data
        raw_df = data.get("raw_data")
        if raw_df is not None and not raw_df.empty:
            # Filter selectors
            # Sex / gender
            genders = sorted(raw_df["sex"].dropna().unique().tolist())
            selected_genders = st.multiselect(
                "Filter by sex", options=genders, default=["M"]
            )
            # Response
            responses = sorted(raw_df["response"].dropna().unique().tolist())
            selected_responses = st.multiselect(
                "Filter by response", options=responses, default=["yes"]
            )
            # Condition
            conditions = sorted(raw_df["condition"].dropna().unique().tolist())
            selected_conditions = st.multiselect(
                "Filter by condition", options=conditions, default=["melanoma"]
            )
            # Time from treatment start
            times = sorted(
                raw_df["time_from_treatment_start"].dropna().unique().tolist()
            )
            selected_times = st.multiselect(
                "Filter by time from treatment start", options=times, default=[0]
            )
            # Apply filters
            filt_df = raw_df[
                (raw_df["sex"].isin(selected_genders))
                & (raw_df["response"].isin(selected_responses))
                & (raw_df["condition"].isin(selected_conditions))
                & (raw_df["time_from_treatment_start"].isin(selected_times))
            ]
            # Display the filtered data
            st.dataframe(
                filt_df.reset_index(drop=True),
                use_container_width=True,
                height=500,
            )
        else:
            st.info("No raw data available to display.")
    elif page == "Responder Analysis":
        st.header("Responder vs Non‑Responder Comparison")
        st.markdown(
            """
            The plots below compare the relative frequencies of immune cell
            populations between responders and non‑responders among
            melanoma patients treated with miraclib.  Use the controls
            to select which populations to display.  Statistical
            significance (Welch’s t‑test) is indicated in the table.
            """
        )
        # Filter by population for plotting
        populations = sorted(stats_df["population"].tolist())
        selected_pop = st.multiselect(
            "Select populations to plot", options=populations, default=populations
        )
        
        # Build a combined DataFrame for plotting
        # We'll query from database: join samples and cell_counts to compute percentages
        db_path_obj = Path(db_path)
        conn = sqlite3.connect(db_path_obj)
        query = """
            SELECT s.sample, s.response, c.population, c.count
            FROM samples s
            JOIN cell_counts c ON s.id = c.sample_id
            WHERE s.condition = 'melanoma'
              AND s.sample_type = 'PBMC'
              AND s.treatment = 'miraclib'
        """
        df = pd.read_sql_query(query, conn)
        conn.close()
        # Compute total counts per sample
        totals = df.groupby("sample")["count"].sum().rename("total_count")
        df = df.join(totals, on="sample")
        df["percentage"] = df["count"] / df["total_count"] * 100
        df = df[df["population"].isin(selected_pop)]
        # Plot: population vs percentage by response
        if not df.empty:
            fig = px.box(
                df,
                x="population",
                y="percentage",
                color="response",
                title="Relative frequencies by response",
                labels={
                    "percentage": "Relative frequency (%)",
                    "population": "Population",
                },
                points="all",
            )
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No data available for the selected populations.")
        # Show statistics table
        st.subheader("Statistical Summary")
        st.write(
            "Mean relative frequencies and p values calculated using Welch’s t‑test.  Significant differences (p < 0.05) are highlighted."
        )
        stats_display = stats_df.copy()
        if selected_pop:
            stats_display = stats_display[
                stats_display["population"].isin(selected_pop)
            ]
        # Format columns
        stats_display = stats_display.round(
            {
                "mean_percentage_responders": 2,
                "mean_percentage_non_responders": 2,
                "difference": 2,
                "p_value": 4,
            }
        )
        # Add significance column as string
        stats_display["significant"] = stats_display["significant"].map(
            {True: "Yes", False: "No"}
        )
        st.dataframe(
            stats_display.set_index("population"),
            use_container_width=True,
            height=400,
        )
    elif page == "Baseline Summary":
        st.header("Baseline Melanoma PBMC Summary")
        st.markdown(
            """
            This view summarises melanoma PBMC samples treated with
            miraclib at baseline (time = 0).  It reports the number of
            samples per project, responder vs non‑responder counts,
            counts by sex, and the mean B cell count for male
            responders at baseline.
            """
        )
        # Display baseline summary
        col1, col2 = st.columns(2)
        with col1:
            st.subheader("Samples per project")
            st.table(
                pd.DataFrame.from_dict(
                    baseline["samples_per_project"], orient="index", columns=["count"]
                )
            )
            st.subheader("Responders vs Non‑responders")
            st.table(
                pd.DataFrame.from_dict(
                    baseline["responders_vs_nonresponders"],
                    orient="index",
                    columns=["count"],
                )
            )
        with col2:
            st.subheader("Sex counts")
            st.table(
                pd.DataFrame.from_dict(
                    baseline["sex_counts"], orient="index", columns=["count"]
                )
            )
            st.subheader("Mean B cell count for male responders (baseline)")
            mean_b = baseline["mean_b_cell_male_responders"]
            if pd.isna(mean_b):
                st.write("No male responders at baseline.")
            else:
                st.write(f"{mean_b:.2f}")


if __name__ == "__main__":
    main()
