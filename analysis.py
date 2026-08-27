import csv
import sqlite3
from pathlib import Path
from statistics import median

from scipy.stats import mannwhitneyu

ROOT_DIR = Path(__file__).resolve().parent
DB_PATH = ROOT_DIR / "cell-count.db"
OUTPUT_DIR = ROOT_DIR / "output"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


def get_connection(db_path: Path = DB_PATH) -> sqlite3.Connection:
    if not db_path.exists():
        raise FileNotFoundError(
            f"Database not found at {db_path}. Run `python load_data.py` first."
        )
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    return conn


def relative_frequencies(conn: sqlite3.Connection) -> list[dict]:
    """Part 2: relative frequency (%) of each cell population within each sample."""
    query = """
        SELECT
            sm.sample_code AS sample,
            SUM(cc.count) OVER (PARTITION BY cc.sample_id) AS total_count,
            cc.population AS population,
            cc.count AS count,
            ROUND(
                100.0 * cc.count / SUM(cc.count) OVER (PARTITION BY cc.sample_id), 2
            ) AS percentage
        FROM cell_counts cc
        JOIN samples sm ON sm.sample_id = cc.sample_id
        ORDER BY sm.sample_code, cc.population
    """
    return [dict(row) for row in conn.execute(query)]


def melanoma_miraclib_pbmc_frequencies(conn: sqlite3.Connection) -> list[dict]:
    """Part 3 (plot data): per-sample relative frequencies (%) for melanoma
    patients on miraclib, PBMC samples only, labeled by responder status.
    One row per (sample, population) -- shape a boxplot grouped by
    `response` and faceted by `population` can consume directly.
    """
    query = """
        SELECT
            sm.sample_code AS sample,
            s.response AS response,
            cc.population AS population,
            ROUND(
                100.0 * cc.count / SUM(cc.count) OVER (PARTITION BY cc.sample_id), 2
            ) AS percentage
        FROM cell_counts cc
        JOIN samples sm ON sm.sample_id = cc.sample_id
        JOIN subjects s ON s.subject_id = sm.subject_id
        WHERE s.condition = 'melanoma'
          AND s.treatment = 'miraclib'
          AND sm.sample_type = 'PBMC'
        ORDER BY cc.population, s.response, sm.sample_code
    """
    return [dict(row) for row in conn.execute(query)]


def _benjamini_hochberg(p_values: list[float]) -> list[float]:
    """FDR-adjusted p-values (Benjamini-Hochberg step-up procedure)."""
    m = len(p_values)
    order = sorted(range(m), key=lambda i: p_values[i])
    adjusted = [0.0] * m
    running_min = 1.0
    for rank in range(m - 1, -1, -1):
        idx = order[rank]
        value = p_values[idx] * m / (rank + 1)
        running_min = min(running_min, value)
        adjusted[idx] = min(running_min, 1.0)
    return adjusted


def _responder_significance_from_frequencies(freq_rows: list[dict]) -> list[dict]:
    """Mann-Whitney U (two-sided) per population, BH-adjusted across populations."""
    by_population: dict[str, dict[str, list[float]]] = {
        pop: {"yes": [], "no": []} for pop in POPULATIONS
    }
    for row in freq_rows:
        by_population[row["population"]][row["response"]].append(row["percentage"])

    results = []
    for population in POPULATIONS:
        responders = by_population[population]["yes"]
        non_responders = by_population[population]["no"]
        stat, p_value = mannwhitneyu(responders, non_responders, alternative="two-sided")
        results.append(
            {
                "population": population,
                "n_responders": len(responders),
                "n_non_responders": len(non_responders),
                "median_responders": round(median(responders), 2),
                "median_non_responders": round(median(non_responders), 2),
                "u_statistic": stat,
                "p_value": p_value,
            }
        )

    adjusted = _benjamini_hochberg([r["p_value"] for r in results])
    for row, p_adj in zip(results, adjusted):
        row["p_value"] = round(row["p_value"], 6)
        row["p_value_adj"] = round(p_adj, 6)
        row["significant"] = p_adj < 0.05

    results.sort(key=lambda r: r["p_value_adj"])
    return results


def responder_significance(conn: sqlite3.Connection) -> list[dict]:
    """Part 3b (stats): for melanoma/miraclib/PBMC samples, test whether each
    cell population's relative frequency differs between responders and
    non-responders. NOTE: pools all of a subject's samples (day 0, 7, 14) as if
    they were independent observations, which is not strictly correct -- see
    `baseline_responder_significance` for the repeated-measures-safe version.
    """
    return _responder_significance_from_frequencies(melanoma_miraclib_pbmc_frequencies(conn))


BASELINE_MELANOMA_MIRACLIB_PBMC_FILTER = """
    s.condition = 'melanoma'
    AND s.treatment = 'miraclib'
    AND sm.sample_type = 'PBMC'
    AND sm.time_from_treatment_start = 0
"""


def baseline_responder_frequencies(conn: sqlite3.Connection) -> list[dict]:
    """Part 3c (plot data): same comparison as `melanoma_miraclib_pbmc_frequencies`,
    restricted to baseline (day 0) samples so each subject contributes exactly
    one observation per population -- avoids treating a subject's day 0/7/14
    samples as independent when they are repeated measures on the same person.
    """
    query = f"""
        SELECT
            sm.sample_code AS sample,
            s.response AS response,
            cc.population AS population,
            ROUND(
                100.0 * cc.count / SUM(cc.count) OVER (PARTITION BY cc.sample_id), 2
            ) AS percentage
        FROM cell_counts cc
        JOIN samples sm ON sm.sample_id = cc.sample_id
        JOIN subjects s ON s.subject_id = sm.subject_id
        WHERE {BASELINE_MELANOMA_MIRACLIB_PBMC_FILTER}
        ORDER BY cc.population, s.response, sm.sample_code
    """
    return [dict(row) for row in conn.execute(query)]


def baseline_responder_significance(conn: sqlite3.Connection) -> list[dict]:
    """Part 3c (stats): same test as `responder_significance`, but on baseline-only
    samples -- one independent observation per subject, so the independence
    assumption behind Mann-Whitney U actually holds.
    """
    return _responder_significance_from_frequencies(baseline_responder_frequencies(conn))


def baseline_melanoma_miraclib_pbmc_samples(conn: sqlite3.Connection) -> list[dict]:
    """Part 4: baseline (day 0) melanoma/miraclib/PBMC samples."""
    query = f"""
        SELECT
            p.project_name AS project,
            s.subject_code AS subject,
            s.sex AS sex,
            s.response AS response,
            sm.sample_code AS sample
        FROM samples sm
        JOIN subjects s ON s.subject_id = sm.subject_id
        JOIN projects p ON p.project_id = s.project_id
        WHERE {BASELINE_MELANOMA_MIRACLIB_PBMC_FILTER}
        ORDER BY p.project_name, s.subject_code
    """
    return [dict(row) for row in conn.execute(query)]


def baseline_subset_summary(conn: sqlite3.Connection) -> list[dict]:
    """Part 4: counts within the baseline subset -- samples per project,
    subjects per responder status, subjects per sex.
    """
    by_project = f"""
        SELECT 'project' AS category, p.project_name AS value, COUNT(*) AS count
        FROM samples sm
        JOIN subjects s ON s.subject_id = sm.subject_id
        JOIN projects p ON p.project_id = s.project_id
        WHERE {BASELINE_MELANOMA_MIRACLIB_PBMC_FILTER}
        GROUP BY p.project_name
    """
    by_response = f"""
        SELECT 'response' AS category, s.response AS value, COUNT(DISTINCT s.subject_id) AS count
        FROM samples sm
        JOIN subjects s ON s.subject_id = sm.subject_id
        WHERE {BASELINE_MELANOMA_MIRACLIB_PBMC_FILTER}
        GROUP BY s.response
    """
    by_sex = f"""
        SELECT 'sex' AS category, s.sex AS value, COUNT(DISTINCT s.subject_id) AS count
        FROM samples sm
        JOIN subjects s ON s.subject_id = sm.subject_id
        WHERE {BASELINE_MELANOMA_MIRACLIB_PBMC_FILTER}
        GROUP BY s.sex
    """
    rows = []
    for query in (by_project, by_response, by_sex):
        rows.extend(dict(row) for row in conn.execute(query))
    return rows

def google_form_query(conn: sqlite3.Connection) -> list[dict]:
    query = f"""
        SELECT AVG(cc.count) AS avg_b_cell_count
        FROM cell_counts cc
        JOIN samples sm ON sm.sample_id = cc.sample_id
        JOIN subjects s ON s.subject_id = sm.subject_id
        WHERE s.condition = 'melanoma'
        AND s.sex = 'M'
        AND s.response = 'yes'
        AND sm.time_from_treatment_start = 0
        AND cc.population = 'b_cell'
    """
    return [dict(row) for row in conn.execute(query)]


QUESTIONS = {
    "2": ("Relative frequency of each cell population per sample", relative_frequencies),
    "3a": (
        "Melanoma/miraclib/PBMC relative frequencies by responder status (plot data)",
        melanoma_miraclib_pbmc_frequencies,
    ),
    "3b": (
        "Significance of responder vs non-responder frequency differences per population",
        responder_significance,
    ),
    "3c": (
        "Significance of responder vs non-responder differences, baseline-only "
        "(one independent sample per subject)",
        baseline_responder_significance,
    ),
    "4a": (
        "Baseline melanoma/miraclib/PBMC samples",
        baseline_melanoma_miraclib_pbmc_samples,
    ),
    "4b": (
        "Baseline subset counts by project, response, and sex",
        baseline_subset_summary,
    ),
}


def write_csv(rows: list[dict], path: Path) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    conn = get_connection()
    try:
        for key, (description, func) in QUESTIONS.items():
            rows = func(conn)
            out_path = OUTPUT_DIR / f"part_{key}.csv"
            write_csv(rows, out_path)
            print(f"Part {key}: {description} -> wrote {len(rows)} rows to {out_path}")

        google_rows = google_form_query(conn)
        out_path = OUTPUT_DIR / "google_form_query.csv"
        write_csv(google_rows, out_path)
        print(f"Wrote {len(google_rows)} rows to {out_path}")
    finally:
        conn.close()


if __name__ == "__main__":
    main()
