# Immune Cell Population Analysis

This project consists of two main components:

1. **Data management and analysis scripts** that construct a
   relational database from the provided CSV file, compute summary
   statistics and perform statistical tests.
2. **An interactive dashboard** built with [Streamlit](https://streamlit.io/)
   that allows exploration of the data, visualisation of responder
   differences and inspection of baseline characteristics.

## Repository structure

```
immune-cell-dashboard/
├── data_management.py  # create SQLite database and load CSV
├── analysis.py         # compute summary statistics and run tests
├── app.py              # Streamlit dashboard
├── immune_data.db      # SQLite database created from the CSV (generated at runtime)
├── README.md           # this file
└── requirements.txt    # Python dependencies
```

## Database schema

The database schema is deliberately normalised to accommodate hundreds of
projects, thousands of samples and additional cell types without
restructuring.  It comprises two tables:

### `samples`

| Column | Type | Description |
|-------|------|-------------|
| **id** | INTEGER PRIMARY KEY | Surrogate key assigned when loading data |
| **sample** | TEXT | Sample identifier matching the `sample` column in the CSV |
| **subject** | TEXT | Subject identifier |
| **project** | TEXT | Project name |
| **condition** | TEXT | Disease or indication (e.g. `melanoma`, `carcinoma`) |
| **age** | INTEGER | Age of the subject |
| **sex** | TEXT | Sex of the subject (`M` or `F`) |
| **treatment** | TEXT | Name of treatment administered |
| **response** | TEXT | Response to treatment (`yes` for responders, `no` for non‑responders) |
| **sample_type** | TEXT | Type of sample (e.g. `PBMC`) |
| **time_from_treatment_start** | INTEGER | Time elapsed from treatment start |

### `cell_counts`

| Column | Type | Description |
|-------|------|-------------|
| **id** | INTEGER PRIMARY KEY | Surrogate key |
| **sample_id** | INTEGER | Foreign key referencing `samples.id` |
| **population** | TEXT | Name of the immune cell population (e.g. `b_cell`) |
| **count** | INTEGER | Count for that population |

In this normalised design the `samples` table contains **only the
metadata** for each biological sample and does not store the cell
counts.  All immune cell population counts are stored in the
`cell_counts` table along with a reference to the sample they belong to.
If new immune populations are added, they can be loaded by inserting
additional rows into `cell_counts` without altering the structure of
the `samples` table.  For large datasets, indexes on `project`,
`sample_type`, `population` or `time_from_treatment_start` can
accelerate queries for specific analyses.

### Design rationale

This schema follows the principles of relational database design, in
particular *Third Normal Form* (3NF).  In the original CSV the counts
for each immune population were represented as separate columns on
every row.  That flat, wide layout is convenient for spreadsheets but
problematic for a database: adding a new cell type would require
altering the table structure, and repeating population columns across
rows introduces redundancy and invites update anomalies.  By separating
the cell counts into their own table and storing one row per
(sample, population) pair, the design treats the immune population
names as *data* rather than part of the schema.  This makes the system
extensible (new populations simply add rows), reduces redundancy
(`samples` stores metadata once per sample), and simplifies
analytics.

## Run

1. Install the dependencies (preferably in a virtual environment):

   ```bash
   pip install -r requirements.txt
   ```

2. Initialise the database and load the CSV data.  Assuming the CSV is
   located at `../cell-count.csv` relative to the repository root,
   execute:

   ```bash
   python data_management.py ./cell-count.csv immune_data.db
   ```

   This creates `immune_data.db` in the `immune-cell-dashboard` directory
   containing the normalised tables and all 10 500 samples.

3. To run the command‑line analysis scripts directly (optional):

   ```bash
   python analysis.py immune_data.db
   ```

   This prints the relative frequency table (first few rows), the
   responder analysis and the baseline summary to standard output.

4. To start the interactive dashboard:

   ```bash
   streamlit run app.py
   ```
   - **Data Overview:** Explore relative frequencies for each sample and
     population.  Use filters to focus on particular samples or cell
     types.
   - **Responder Analysis:** Visualise boxplots comparing relative
     frequencies between responders and non‑responders among melanoma
     PBMC samples treated with miraclib.  A table lists mean
     percentages and p values from Welch’s t‑tests with significant
     differences highlighted.
   - **Baseline Summary:** Display counts of baseline melanoma PBMC
     samples (time = 0) per project, responder status and sex.  It also
     reports the mean number of B cells for male responders at baseline.

5. You can also create a GitHub Codespace to run:
    ```bash
   pip install -r immune-cell-dashboard/requirements.txt
   python immune-cell-dashboard/data_management.py ./cell-count.csv immune_data.db
   streamlit run immune-cell-dashboard/app.py --server.address 0.0.0.0 --server.port 8501
   ```

## Code structure rationale

The code is organised into separate modules to encourage clarity and
reusability:

* **data_management.py** contains everything related to the database:
  schema definition and a loader that ingests the original CSV.  The
  functions are idempotent and can be called independently.
* **analysis.py** implements the core analytics: computing relative
  frequencies, performing responder comparisons and calculating
  baseline summaries.  These functions operate purely on the SQLite
  database and return pandas DataFrames or Python dicts, making them
  amenable to both command‑line use and interactive consumption.
* **app.py** stitches everything together into a Streamlit dashboard.
  It caches the heavy database queries, provides an intuitive user
  interface for exploring the data, and displays plots and tables.

This modularity means that if Bob wants to add more analytics later
(e.g. new statistical tests or different filtering criteria) he can
implement them in `analysis.py` and immediately surface them in the
Streamlit app.
