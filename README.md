# UrbanScope

**UrbanScope** is a fully automated, GitHub-native system for discovering, tracking, and publishing
**urban metagenomics and metatranscriptomics studies and their associated datasets**.

It continuously scans public biomedical repositories for new urban omics studies, links publications
to sequencing datasets (PubMed ↔ SRA), enriches records with location and study-type metadata,
and publishes a live, searchable dashboard and downloadable archives — all at **zero cost**.

---

## 🌍 Live Dashboard

```
https://aglucaci.github.io/urbanscope/
```

Includes:
- latest newly discovered datasets
- interactive filtering (study type, city, country)
- archive browser by year
- CSV downloads for collaborators

---

## 🎯 What UrbanScope Does

UrbanScope automatically:

- Detects **urban / built-environment / wastewater / air / surface** studies
- Identifies **metagenomics & metatranscriptomics** datasets
- Links:
  - **PubMed papers**
  - **SRA datasets**
- Enriches records with:
  - **Study type**: `air`, `wastewater`, `surface`, `other`
  - **City / country** (best-effort heuristic extraction)
- Maintains a **permanent, deduplicated catalog**
- Publishes:
  - JSON (machine-readable)
  - CSV (collaborator-ready)
  - HTML dashboard (GitHub Pages)

No APIs, no databases, no cloud accounts.

---

## 🧬 Data Sources

- **NCBI PubMed** (literature discovery)
- **NCBI SRA** (sequencing datasets)
- Linking via **NCBI E-utilities (elink)**

All data are public and accessed using NCBI-approved rate limits.

---

## 🧠 Design Philosophy

- **Situational awareness, not retrospection**
- **Append-only, auditable records**
- **Static artifacts over fragile backends**
- **GitHub as compute + storage + publishing**
- **Free forever**

UrbanScope is intended as **scientific infrastructure**, not a demo.

---

## 📁 Repository Structure

```
urbanscope/
├── scripts/
│   └── urbanscope_radar.py
├── data/
│   ├── seen_ids.txt
│   ├── catalog_2016.jsonl
│   ├── catalog_2017.jsonl
│   └── ...
├── docs/
│   ├── index.html
│   ├── latest.json
│   └── archive/
│       ├── index.html
│       ├── index.json
│       └── csv/
│           ├── catalog_2016.csv
│           ├── catalog_2017.csv
│           └── latest_added.csv
└── .github/
    └── workflows/
        ├── backfill-year.yml
        └── daily.yml
```

---

## ⚙️ How It Runs (GitHub-Only)

### Historical Backfill (Year-by-Year)
Run manually from **GitHub Actions** using the workflow:

```
UrbanScope Backfill (Year)
```

Each run:
- Appends to `data/catalog_<YEAR>.jsonl`
- Updates CSV exports
- Updates the public archive

---

### Daily Incremental Updates
Runs automatically every day via cron:

- Detects **new studies only**
- Deduplicates permanently
- Appends to the current year’s catalog
- Updates dashboard and CSVs

---

## 📊 Output Formats

- **JSONL** — append-only archival storage
- **JSON** — live dashboard updates
- **CSV** — collaborator-friendly downloads
- **HTML** — static dashboard (GitHub Pages)

All outputs are versioned in Git.

---

## 🔒 Cost & Infrastructure

| Component | Provider | Cost |
|--------|---------|------|
| Compute | GitHub Actions | Free |
| Storage | GitHub repo | Free |
| Hosting | GitHub Pages | Free |
| APIs | NCBI E-utilities | Free |

**No credit card required.**

---

## ⚠️ Disclaimer

UrbanScope is provided for **research and informational purposes only**.
It does not constitute public-health guidance, policy recommendations,
or medical advice.

---

## 👤 Author

**Alexander G. Lucaci, PhD**  
Computational Evolutionary Biology • Urban Metagenomics • Genomic Surveillance
