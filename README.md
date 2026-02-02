# UrbanScope SRR Harvester (Modular)

This is a modular refactor of the SRR harvester with **hard output size limits**:
- Any JSONL output is **rotated** so each file stays <= **50MB**
- Large JSON exports are written as **chunked JSON array parts** + a manifest

## Run

From repo root:

```bash
python -m scripts.urbanscope_harvester daily --recent-days 7 --max-per-day 500 --fetch-biosample --fetch-bioproject
python -m scripts.urbanscope_harvester crawl --page-size 500 --sort date --fetch-biosample --fetch-bioproject
python -m scripts.urbanscope_harvester backfill-year --year 2024 --max-per-day 500

# What I ran 

python -m scripts.urbanscope_harvester crawl --page-size 100 --sort date --fetch-biosample --fetch-bioproject --max-total 10000

```

## Outputs

- `data/srr_catalog_<YEAR>.jsonl` (+ rotated parts `_partNNN.jsonl`)
- `docs/db/srr_records_partNNN.json` + `docs/db/srr_records_manifest.json`
- `docs/latest_srr.json` (trimmed if needed to stay <= 50MB)

## Notes
- Set `NCBI_API_KEY` env var to be nice to NCBI and reduce throttling.
- Set `NCBI_EMAIL` for proper User-Agent identification.



python -m http.server 8000
http://localhost:8000