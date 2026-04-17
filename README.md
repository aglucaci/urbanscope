# UMDB (Urban Microbiome Database)

**UMDB (Urban Microbiome Database)** is a GitHub-hosted bioinformatics resource for discovering, organizing, and publishing urban environmental sequencing metadata from public repositories.

The repository has two jobs:

1. Harvest and normalize public metadata from NCBI-linked resources.
2. Publish those records as a static website and downloadable database under `docs/`.

The current public site is available at:

https://aglucaci.github.io/umdb/

## What This Repository Is

UMDB (Urban Microbiome Database) is designed as scientific infrastructure for urban metagenomics and related environmental sequencing studies. It focuses on public sequencing records associated with urban, built-environment, wastewater, air, surface, and similar city-linked contexts.

Rather than operating as a traditional server-backed application, the repository stores harvested records directly in versioned files and publishes a static web interface through GitHub Pages. This keeps the project easy to inspect, easy to archive, and inexpensive to maintain.

## What The Website Provides

The published site includes:

- a searchable BioProject explorer
- expandable SRR-level detail views
- filter controls for geography, assay class, center, and year
- search-derived dataset bundling for matched SRR cohorts
- exportable SRR accession lists and FASTQ download scripts
- downloadable JSON artifacts under `docs/db/`
- reviewer-facing pages such as `About`, `Methodology`, and `Data Access`
- summary analytics across countries, cities, years, centers, and assay classes
- a global map view of country-level coverage

## How The Repository Works

At a high level, UMDB (Urban Microbiome Database) follows this flow:

1. Discover candidate public records from NCBI-linked resources.
2. Retrieve run-level metadata for SRR entries.
3. Join associated BioProject and BioSample context when available.
4. Derive practical annotations such as assay class and geographic labels.
5. Store the resulting records in local append-only or chunked artifacts.
6. Export static JSON files that power the website.
7. Publish the `docs/` directory through GitHub Pages.

This means the repository itself is both the processing workspace and the published database release.

## Data Sources

UMDB (Urban Microbiome Database) uses public repository metadata. The main sources currently represented in the code and outputs are:

- **NCBI SRA** for run-level sequencing metadata
- **NCBI BioProject** for project-level metadata
- **NCBI BioSample** for sample and location context
- **NCBI E-utilities** for public metadata retrieval and cross-linking

All records remain subject to the quality and completeness of the original public submissions.

## Discovery Strategy

UMDB (Urban Microbiome Database) no longer relies on a single search string for study discovery. The harvester now supports a multi-profile query strategy intended to improve recall across different ways urban microbiome studies are described in public repositories.

The default discovery pass combines several query profiles:

- urban shotgun and metagenomic studies
- urban amplicon and marker-gene studies
- wastewater and sewage surveillance studies
- transit and surface microbiome studies
- urban air, aerosol, and dust studies

This improves coverage for studies that do not all use the same vocabulary in titles, abstracts, or repository metadata.

Even with broader discovery, UrbanScope should still be treated as a high-recall public index rather than a guaranteed complete census of every urban microbiome study ever performed. Repository metadata are inconsistent, and some studies are only discoverable through follow-up curation.

## AI-Assisted Curation

UrbanScope can now add an AI-assisted curation layer on top of the original public metadata.

This curation layer is designed to:

- judge whether the available metadata are sufficient for city-level interpretation
- decide whether the sample appears to come from a truly urban or city-associated context
- inspect BioSample attributes, titles, descriptions, and run metadata together
- propose corrected country, city, or assay annotations when the evidence supports it
- preserve the original annotations while publishing AI-reviewed final values alongside them

AI curation is stored separately in cache and merged into the published dataset at export time as an `ai_curation` object on each record. The original repository metadata remain unchanged.

### AI Curation Fields

When present, `ai_curation` can include fields such as:

- `metadata_sufficiency`
- `metadata_sufficient`
- `urban_origin`
- `confidence`
- `environment_type`
- `ai_fixed`
- `original_country`
- `original_city`
- `original_assay_class`
- `final_country`
- `final_city`
- `final_assay_class`
- `reasoning_summary`
- `evidence`

The website can then distinguish between original-only records and AI-fixed records.

## Repository Layout

The repository is organized around three main areas: harvesting code, intermediate data, and published web artifacts.

```text
urbanscope/
├── data/
│   ├── cache/
│   │   ├── bioproject.json
│   │   ├── bioproject_uid.json
│   │   └── biosample.json
│   ├── seen_sra_uids.txt
│   ├── seen_srr_runs.txt
│   └── srr_catalog_2026*.jsonl
├── docs/
│   ├── index.html
│   ├── about.html
│   ├── methodology.html
│   ├── data.html
│   ├── analytics.html
│   ├── global-map.html
│   ├── latest_srr.json
│   ├── assets/
│   │   ├── db.js
│   │   └── site.css
│   └── db/
│       ├── srr_records_manifest.json
│       ├── srr_records_part000.json
│       ├── ...
│       ├── srr_index.json
│       ├── bioprojects.json
│       └── biosamples.json
├── logo/
├── scripts/
│   └── urbanscope_harvester/
│       ├── cli.py
│       ├── ingest.py
│       ├── exports.py
│       ├── ncbi.py
│       ├── bioproject.py
│       ├── biosample.py
│       ├── assay.py
│       ├── utils.py
│       └── config.py
└── README.md
```

## What Lives In `scripts/urbanscope_harvester/`

The Python package under `scripts/urbanscope_harvester/` is the operational core of the project.

- `cli.py` defines command-line entry behavior.
- `ingest.py` coordinates harvesting and record assembly.
- `ncbi.py` contains public NCBI request logic.
- `bioproject.py` and `biosample.py` handle project and sample enrichment.
- `assay.py` derives higher-level assay labels from run metadata.
- `exports.py` writes the static JSON artifacts used by the website.
- `config.py` centralizes path and output settings.
- `utils.py` provides common helpers for JSON, timestamps, and iteration.

Together, these modules take public metadata and turn it into the chunked release files under `docs/db/`.

## What Lives In `data/`

The `data/` directory is the working data store for harvesting and enrichment.

- `seen_sra_uids.txt` and `seen_srr_runs.txt` track identifiers that have already been processed.
- `srr_catalog_*.jsonl` files store accumulated run-level records.
- `cache/` stores locally reused metadata lookups such as BioProject and BioSample responses.

This directory is important because it represents the stateful layer of the pipeline. It is where UrbanScope remembers what it has seen and where it stages enriched metadata before export.

## What Lives In `docs/`

The `docs/` directory is the public release of the project.

- HTML pages provide the user-facing web interface.
- `docs/db/` contains machine-readable database artifacts.
- `latest_srr.json` provides a compact recent-record feed.
- `assets/` contains shared JavaScript and CSS used across the site.

Anything under `docs/` should be treated as published output, not just internal project files.

## Database Model

UrbanScope currently publishes a manifest-plus-parts layout for SRR-level records.

### Manifest

`docs/db/srr_records_manifest.json` describes:

- when the release was generated
- how many total records exist
- which part files belong to the release
- which years are represented in the export summary

### Part Files

`docs/db/srr_records_partXXX.json` files contain arrays of SRR-level records. A typical record can include:

- `srr` for the run accession
- `runinfo_row` for SRA RunInfo fields
- `bioproject` for project-level metadata
- `geo` for parsed geographic context
- `assay` for high-level assay classification
- `ncbi` for source URLs

The website loads these files directly in the browser and aggregates them into BioProject-level summaries for exploration.

## Why The Site Aggregates By BioProject

Many users want to discover studies, not just individual runs. UrbanScope therefore groups records by BioProject in the main explorer while preserving the underlying SRR records for inspection.

This gives the site two layers:

- a study-facing summary layer for browsing
- a run-facing provenance layer for detail and reuse

That structure is useful for both casual exploration and manuscript review.

## Publication And Reviewer-Facing Pages

The site now includes pages that help present the repository as a maintained scientific database:

- `about.html` explains the scope and intended use of the resource
- `methodology.html` describes derivation logic, provenance, and caveats
- `data.html` explains downloads and schema expectations
- `analytics.html` summarizes global counts and distributions
- `global-map.html` visualizes country-level coverage

These pages are meant to reduce ambiguity for collaborators, manuscript reviewers, and future users of the database.

## Design Principles

UrbanScope follows a few simple design principles:

- static publication over server complexity
- auditable files over opaque backend state
- public metadata reuse over proprietary infrastructure
- browser access plus machine-readable exports
- low-cost maintenance with transparent artifacts

## Typical Update Pattern

Although the exact command flow may evolve, a normal update cycle looks like this:

1. Harvest new or target records into `data/`.
2. Update caches and enrichments.
3. Rebuild exported JSON release artifacts in `docs/db/`.
4. Refresh site pages that consume the latest summaries.
5. Commit the updated data and published files.

In other words, a repository update is also a database release.

## Local AI Curation

To use AI-assisted metadata curation locally, set your environment first:

```bash
export OPENAI_API_KEY="your_openai_api_key"
export OPENAI_MODEL="gpt-4o-mini"
export NCBI_API_KEY="your_ncbi_api_key"
export NCBI_EMAIL="you@example.org"
export NCBI_TOOL="urbanscope-srr-harvester"
```

You can curate newly harvested records during a deep crawl:

```bash
python3 -m scripts.urbanscope_harvester.cli crawl \
  --fetch-biosample \
  --fetch-bioproject \
  --page-size 500 \
  --sort date \
  --max-total 0 \
  --stop-after-new-srr 0 \
  --ai-curate
```

You can also run AI curation against existing harvested records:

```bash
python3 -m scripts.urbanscope_harvester.cli curate-ai \
  --model gpt-4o-mini \
  --max-records 500
```

To overwrite previously cached AI reviews:

```bash
python3 -m scripts.urbanscope_harvester.cli curate-ai \
  --model gpt-4o-mini \
  --overwrite \
  --max-records 500
```

To target a single year:

```bash
python3 -m scripts.urbanscope_harvester.cli curate-ai \
  --year 2024 \
  --max-records 500
```

## GitHub Actions AI Curation

The repository also includes a dedicated full-dataset AI curation workflow:

- Workflow file: `.github/workflows/urbanscope_full_ai_curation.yml`
- Workflow name in GitHub Actions: `UrbanScope Full AI Curation`
- Required secret: `OPENAI_API_KEY`
- Optional but recommended secrets for NCBI-linked context: `NCBI_API_KEY`, `NCBI_EMAIL`

To run it in GitHub:

1. Open the repository Actions tab.
2. Select `UrbanScope Full AI Curation`.
3. Click `Run workflow`.
4. Leave `max_records` as `0` to process the entire dataset.
5. Leave `year` blank to process all years, or set a single year for a targeted pass.
6. Turn on `overwrite` only when you want to replace cached AI reviews.

This workflow runs `python3 -m scripts.urbanscope_harvester.cli curate-ai`, rebuilds the exported database artifacts, and commits updated `data/` and `docs/` outputs back to the repository.

## Output Formats

UrbanScope publishes several output styles because different users need different levels of structure:

- **JSONL** for accumulation-oriented internal catalogs in `data/`
- **JSON** for the public website and machine-readable exports in `docs/`
- **HTML** for the user-facing static interface
- **shell scripts and accession lists** generated from browser searches for raw FASTQ retrieval

At the moment, the public website is centered on JSON-driven delivery rather than a live backend database.

## Downloading Raw Datasets

UrbanScope is primarily a metadata and discovery resource, but the web explorer now supports bundling matched search results into download-ready dataset cohorts.

From the main explorer, a user can:

- export a JSON bundle of the currently matched runs
- export a plain-text list of SRR accessions
- export a shell script that uses `prefetch` and `fasterq-dump` from SRA Toolkit

This makes it possible to go from a filtered UrbanScope search to raw FASTQ retrieval without manually collecting accessions one by one.

## Limitations

UrbanScope is useful, but it should be interpreted carefully.

- Geographic metadata may be incomplete, ambiguous, or inconsistently formatted.
- Assay labels are practical summary categories, not formal reannotation of every study.
- Inclusion reflects what is publicly available and discoverable through the current harvesting logic.
- Presence in UrbanScope does not imply endorsement, quality ranking, or uniform study design.

These limitations are expected for public metadata integration projects and should be stated clearly in any manuscript built around the resource.

## Cost And Infrastructure Model

The project is intentionally lightweight.

| Component | Provider | Cost |
| --- | --- | --- |
| Compute | Local runs or GitHub-based workflows | Low to free |
| Storage | Git repository artifacts | Low to free |
| Hosting | GitHub Pages | Free |
| Metadata access | NCBI public utilities | Free within usage guidance |

## Disclaimer

UrbanScope is provided for research and informational purposes only. It does not constitute medical advice, clinical guidance, or public health policy.

## Author

**Alexander G. Lucaci, PhD**  
Computational Evolutionary Biology  
Urban Metagenomics  
Genomic Surveillance
