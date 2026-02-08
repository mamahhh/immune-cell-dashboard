from __future__ import annotations

import csv
import sqlite3
from pathlib import Path


def initialize_database(db_path: str | Path) -> None:
    """Initialize a SQLite database with the required schema.

    Parameters
    ----------
    db_path : str or Path
        Path to the SQLite database file.  If the file does not exist it
        will be created; if it does exist its schema will be dropped and
        recreated.
    """
    db_path = Path(db_path)
    if db_path.exists():
        # Remove the existing database to start fresh
        db_path.unlink()
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    # Enable foreign key constraints
    c.execute("PRAGMA foreign_keys = ON;")
    # Create tables
    c.execute(
        """
        CREATE TABLE samples (
            id INTEGER PRIMARY KEY,
            sample TEXT NOT NULL,
            subject TEXT NOT NULL,
            project TEXT NOT NULL,
            condition TEXT NOT NULL,
            age INTEGER,
            sex TEXT,
            treatment TEXT,
            response TEXT,
            sample_type TEXT,
            time_from_treatment_start INTEGER
        );
        """
    )
    c.execute(
        """
        CREATE TABLE cell_counts (
            id INTEGER PRIMARY KEY,
            sample_id INTEGER NOT NULL,
            population TEXT NOT NULL,
            count INTEGER NOT NULL,
            FOREIGN KEY(sample_id) REFERENCES samples(id) ON DELETE CASCADE
        );
        """
    )
    conn.commit()
    conn.close()


def load_data(csv_path: str | Path, db_path: str | Path) -> None:
    """Load data from a CSV file into the SQLite database.

    Parameters
    ----------
    csv_path : str or Path
        Path to the CSV file containing cell counts and sample metadata.
    db_path : str or Path
        Path to the SQLite database created by ``initialize_database``.

    The function reads the CSV row by row, inserting each sample's
    metadata into the ``samples`` table and the associated cell counts
    into ``cell_counts``.  It assumes that each row includes the
    columns: ``project``, ``subject``, ``condition``, ``age``,
    ``sex``, ``treatment``, ``response``, ``sample``, ``sample_type``,
    ``time_from_treatment_start``, and the immune cell count columns
    ``b_cell``, ``cd8_t_cell``, ``cd4_t_cell``, ``nk_cell``, and
    ``monocyte``.  The cell count columns are not stored in
    ``samples`` but are normalised into ``cell_counts``.
    """
    csv_path = Path(csv_path)
    db_path = Path(db_path)
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    # Prepare insert statements
    sample_insert = (
        """
        INSERT INTO samples (
            sample, subject, project, condition, age, sex, treatment,
            response, sample_type, time_from_treatment_start
        ) VALUES (?,?,?,?,?,?,?,?,?,?);
        """
    )
    cell_insert = (
        "INSERT INTO cell_counts (sample_id, population, count) VALUES (?,?,?);"
    )
    with csv_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Insert into samples
            sample_values = (
                row["sample"],
                row["subject"],
                row["project"],
                row["condition"],
                int(row["age"]) if row["age"] else None,
                row.get("sex"),
                row.get("treatment"),
                row.get("response"),
                row.get("sample_type"),
                int(float(row.get("time_from_treatment_start", 0))),
            )
            c.execute(sample_insert, sample_values)
            sample_id = c.lastrowid
            # Insert into cell_counts for each population
            populations = [
                ("b_cell", int(row["b_cell"])),
                ("cd8_t_cell", int(row["cd8_t_cell"])),
                ("cd4_t_cell", int(row["cd4_t_cell"])),
                ("nk_cell", int(row["nk_cell"])),
                ("monocyte", int(row["monocyte"])),
            ]
            for population_name, count in populations:
                c.execute(cell_insert, (sample_id, population_name, count))
    conn.commit()
    conn.close()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Initialize a database and load data from a cell-count CSV."
    )
    parser.add_argument(
        "csv_path", type=str, help="Path to the cell-count CSV file"
    )
    parser.add_argument(
        "db_path", type=str, help="Path to the SQLite database to create"
    )
    args = parser.parse_args()
    initialize_database(args.db_path)
    load_data(args.csv_path, args.db_path)