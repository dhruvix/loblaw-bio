import csv
import sqlite3
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parent
CSV_PATH = ROOT_DIR / "cell-count.csv"
DB_PATH = ROOT_DIR / "cell-count.db"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

SCHEMA = """
CREATE TABLE projects (
    project_id   INTEGER PRIMARY KEY AUTOINCREMENT,
    project_name TEXT NOT NULL UNIQUE
);

CREATE TABLE subjects (
    subject_id   INTEGER PRIMARY KEY AUTOINCREMENT,
    subject_code TEXT NOT NULL UNIQUE,
    project_id   INTEGER NOT NULL REFERENCES projects(project_id),
    condition    TEXT,
    age          INTEGER,
    sex          TEXT,
    treatment    TEXT,
    response     TEXT
);

CREATE TABLE samples (
    sample_id                 INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_code               TEXT NOT NULL UNIQUE,
    subject_id                INTEGER NOT NULL REFERENCES subjects(subject_id),
    sample_type               TEXT,
    time_from_treatment_start INTEGER
);

CREATE TABLE cell_counts (
    cell_count_id INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id     INTEGER NOT NULL REFERENCES samples(sample_id),
    population    TEXT NOT NULL CHECK (population IN (
                        'b_cell', 'cd8_t_cell', 'cd4_t_cell', 'nk_cell', 'monocyte'
                   )),
    count         INTEGER NOT NULL,
    UNIQUE (sample_id, population)
);

CREATE INDEX idx_subjects_project ON subjects(project_id);
CREATE INDEX idx_samples_subject ON samples(subject_id);
CREATE INDEX idx_cell_counts_sample ON cell_counts(sample_id);
CREATE INDEX idx_cell_counts_population ON cell_counts(population);
"""


def init_db(conn: sqlite3.Connection) -> None:
    conn.executescript(SCHEMA)


def load_csv(conn: sqlite3.Connection) -> None:
    cur = conn.cursor()
    project_cache: dict[str, int] = {}
    subject_cache: dict[str, int] = {}

    with CSV_PATH.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            project_name = row["project"]
            project_id = project_cache.get(project_name)
            if project_id is None:
                cur.execute(
                    "INSERT INTO projects (project_name) VALUES (?)",
                    (project_name,),
                )
                project_id = cur.lastrowid
                project_cache[project_name] = project_id

            subject_code = row["subject"]
            subject_id = subject_cache.get(subject_code)
            if subject_id is None:
                cur.execute(
                    """INSERT INTO subjects
                       (subject_code, project_id, condition, age, sex, treatment, response)
                       VALUES (?, ?, ?, ?, ?, ?, ?)""",
                    (
                        subject_code,
                        project_id,
                        row["condition"],
                        int(row["age"]),
                        row["sex"],
                        row["treatment"],
                        row["response"],
                    ),
                )
                subject_id = cur.lastrowid
                subject_cache[subject_code] = subject_id

            cur.execute(
                """INSERT INTO samples
                   (sample_code, subject_id, sample_type, time_from_treatment_start)
                   VALUES (?, ?, ?, ?)""",
                (
                    row["sample"],
                    subject_id,
                    row["sample_type"],
                    int(row["time_from_treatment_start"]),
                ),
            )
            sample_id = cur.lastrowid

            cur.executemany(
                "INSERT INTO cell_counts (sample_id, population, count) VALUES (?, ?, ?)",
                [(sample_id, population, int(row[population])) for population in POPULATIONS],
            )

    conn.commit()


def main() -> None:
    if not CSV_PATH.exists():
        raise FileNotFoundError(f"CSV file not found: {CSV_PATH}")

    if DB_PATH.exists():
        DB_PATH.unlink()

    conn = sqlite3.connect(DB_PATH)
    try:
        conn.execute("PRAGMA foreign_keys = ON")
        init_db(conn)
        load_csv(conn)

        counts = {
            table: conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
            for table in ("projects", "subjects", "samples", "cell_counts")
        }
        print(f"Loaded database: {DB_PATH}")
        for table, count in counts.items():
            print(f"  {table}: {count} rows")
    finally:
        conn.close()


if __name__ == "__main__":
    main()
