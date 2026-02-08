from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Dict, Any

import pandas as pd
from scipy import stats


def compute_relative_frequencies(db_path: str | Path) -> pd.DataFrame:
    """Compute the relative frequency of each immune population for each sample.

    Parameters
    ----------
    db_path : str or Path
        Path to the SQLite database created via ``data_management.initialize_database``.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns ``sample``, ``total_count``, ``population``,
        ``count``, and ``percentage`` (relative frequency in percent).  Each
        row corresponds to one sample–population pair.
    """
    db_path = Path(db_path)
    conn = sqlite3.connect(db_path)
    # Query to get counts by sample and population
    df_counts = pd.read_sql_query(
        "SELECT s.sample, c.population, c.count FROM cell_counts c JOIN samples s ON c.sample_id = s.id",
        conn,
    )
    conn.close()
    # Compute total counts per sample
    totals = df_counts.groupby("sample")["count"].sum().rename("total_count")
    df = df_counts.join(totals, on="sample")
    df["percentage"] = (df["count"] / df["total_count"] * 100).astype(float).round(2)
    # Reorder columns
    df = df[["sample", "total_count", "population", "count", "percentage"]]
    return df


def compare_responders(db_path: str | Path) -> pd.DataFrame:
    """Perform statistical comparison between responders and non‑responders.

    Only includes melanoma PBMC samples treated with miraclib.  For
    each immune population, the mean relative frequency is computed
    separately for responders (response = 'yes') and non‑responders
    ('no'), and a two‑sample t‑test (Welch’s t test) is performed on
    the relative frequencies.  The resulting DataFrame contains
    population names, mean percentages for responders and
    non‑responders, the difference of means, t statistics, p values and
    a boolean flag indicating whether the difference is significant at
    α = 0.05.

    Parameters
    ----------
    db_path : str or Path
        Path to the SQLite database.

    Returns
    -------
    pandas.DataFrame
        Summary of the statistical comparison by population.
    """
    db_path = Path(db_path)
    conn = sqlite3.connect(db_path)
    # Query data for melanoma PBMC samples treated with miraclib
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
    df["percentage"] = (df["count"] / df["total_count"] * 100).astype(float).round(2)
    # Group by population and response
    results = []
    for population, group in df.groupby("population"):
        # Split by response
        responders = group[group["response"] == "yes"]["percentage"].values
        non_responders = group[group["response"] == "no"]["percentage"].values
        # Compute means
        mean_resp = (
            responders.mean().astype(float).round(2)
            if len(responders) > 0
            else float("nan")
        )
        mean_non_resp = (
            non_responders.mean().astype(float).round(2)
            if len(non_responders) > 0
            else float("nan")
        )
        # Perform Welch’s t test
        # Use equal_var=False for Welch test; handle potential errors if one group is empty
        if len(responders) > 1 and len(non_responders) > 1:
            t_stat, p_val = stats.ttest_ind(responders, non_responders, equal_var=False)
        else:
            t_stat, p_val = float("nan"), float("nan")
        results.append(
            {
                "population": population,
                "mean_percentage_responders": mean_resp,
                "mean_percentage_non_responders": mean_non_resp,
                "difference": (mean_resp - mean_non_resp).round(2),
                "t_statistic": t_stat.round(4),
                "p_value": p_val.round(4),
                "significant": p_val < 0.05 if not pd.isnull(p_val) else False,
            }
        )
    result_df = pd.DataFrame(results)
    return result_df


def baseline_summary(db_path: str | Path) -> Dict[str, Any]:
    """Summarize baseline melanoma PBMC samples treated with miraclib.

    This function filters the ``samples`` table to include only
    melanoma PBMC samples with treatment ``miraclib`` at
    ``time_from_treatment_start = 0``.  It returns:

    - counts of samples per project
    - counts of responders vs non‑responders
    - counts by sex (M vs F)
    - the mean B cell count for male responders at baseline

    Parameters
    ----------
    db_path : str or Path
        Path to the SQLite database.

    Returns
    -------
    dict
        A dictionary with keys ``samples_per_project``,
        ``responders_vs_nonresponders``, ``sex_counts`` and
        ``mean_b_cell_male_responders``.
    """
    db_path = Path(db_path)
    conn = sqlite3.connect(db_path)
    # Query baseline melanoma PBMC miraclib samples
    # Retrieve baseline melanoma PBMC samples (metadata only)
    baseline_df = pd.read_sql_query(
        """
        SELECT s.id, s.sample, s.project, s.subject, s.condition, s.age, s.sex,
               s.treatment, s.response, s.sample_type, s.time_from_treatment_start
        FROM samples s
        WHERE s.condition = 'melanoma'
          AND s.sample_type = 'PBMC'
          AND s.treatment = 'miraclib'
          AND s.time_from_treatment_start = 0
        """,
        conn,
    )
    # Count by project
    samples_per_project = baseline_df.groupby("project")["sample"].count().to_dict()
    # Count responders vs non responders
    responders_vs_nonresponders = (
        baseline_df.groupby("response")["sample"].count().to_dict()
    )
    # Count by sex
    sex_counts = baseline_df.groupby("sex")["sample"].count().to_dict()
    # Join with cell_counts to get B cell counts for each sample
    b_counts = pd.read_sql_query(
        """
        SELECT sample_id, count AS b_cell_count
        FROM cell_counts
        WHERE population = 'b_cell'
        """,
        conn,
    )
    # Merge baseline metadata with B cell counts
    baseline_df = baseline_df.merge(
        b_counts, left_on="id", right_on="sample_id", how="left"
    )
    # Mean b_cell count for male responders
    male_resp = baseline_df[
        (baseline_df["sex"] == "M") & (baseline_df["response"] == "yes")
    ]
    mean_b = male_resp["b_cell_count"].mean() if not male_resp.empty else float("nan")
    conn.close()
    return {
        "samples_per_project": samples_per_project,
        "responders_vs_nonresponders": responders_vs_nonresponders,
        "sex_counts": sex_counts,
        "mean_b_cell_male_responders": mean_b,
    }


if __name__ == "__main__":
    # Simple demonstration for manual testing
    import argparse

    parser = argparse.ArgumentParser(description="Run analysis on the SQLite database.")
    parser.add_argument("db_path", type=str, help="Path to the SQLite database")
    args = parser.parse_args()
    df_freq = compute_relative_frequencies(args.db_path)
    print("Relative frequencies head:")
    print(df_freq.head())
    df_stats = compare_responders(args.db_path)
    print("Statistical comparison:")
    print(df_stats)
    baseline_info = baseline_summary(args.db_path)
    print("Baseline summary:")
    print(baseline_info)
