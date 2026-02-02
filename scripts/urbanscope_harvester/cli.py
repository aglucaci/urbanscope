from __future__ import annotations
import argparse, datetime as dt
from typing import Any, Dict, List

from .config import (
    DEFAULT_QUERY, BIOSAMPLE_CACHE, BIOPROJECT_CACHE, BIOPROJECT_UID_CACHE,
    DOCS_LATEST_DEBUG, DATA_DIR
)
from .utils import ensure_dirs, load_set, read_json, write_json, append_jsonl
from .ncbi import esearch_recent, esearch_day, esearch_history, esummary_sra
from .ingest import ingest_uids_to_srr, debug_paths
from .exports import rebuild_srr_exports_chunked, write_latest_srr_safe

def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("daily")
    d.add_argument("--days", type=int, default=1)  # kept for compatibility
    d.add_argument("--query", default=DEFAULT_QUERY)
    d.add_argument("--max-per-day", type=int, default=500)
    d.add_argument("--debug", action="store_true")
    d.add_argument("--fetch-biosample", action="store_true")
    d.add_argument("--fetch-bioproject", action="store_true")
    d.add_argument("--runinfo-max-rows", type=int, default=200000)
    d.add_argument("--recent-days", type=int, default=7)

    b = sub.add_parser("backfill-year")
    b.add_argument("--year", type=int, required=True)
    b.add_argument("--query", default=DEFAULT_QUERY)
    b.add_argument("--max-per-day", type=int, default=500)
    b.add_argument("--debug", action="store_true")
    b.add_argument("--fetch-biosample", action="store_true")
    b.add_argument("--fetch-bioproject", action="store_true")
    b.add_argument("--runinfo-max-rows", type=int, default=200000)

    c = sub.add_parser("crawl")
    c.add_argument("--query", default=DEFAULT_QUERY)
    c.add_argument("--page-size", type=int, default=500)
    c.add_argument("--max-total", type=int, default=0)
    c.add_argument("--debug", action="store_true")
    c.add_argument("--fetch-biosample", action="store_true")
    c.add_argument("--fetch-bioproject", action="store_true")
    c.add_argument("--runinfo-max-rows", type=int, default=200000)
    c.add_argument("--stop-after-new-srr", type=int, default=0)
    c.add_argument("--sort", default="date")
    return ap

def run():
    args = build_argparser().parse_args()
    ensure_dirs()

    from .config import SEEN_SRA_UIDS, SEEN_SRR_RUNS, DOCS_LATEST_SRR

    seen_sra = load_set(SEEN_SRA_UIDS)
    seen_srr = load_set(SEEN_SRR_RUNS)

    biosample_cache = read_json(BIOSAMPLE_CACHE, {})
    bp_cache = read_json(BIOPROJECT_CACHE, {})
    bp_uid_cache = read_json(BIOPROJECT_UID_CACHE, {})

    latest_added: List[Dict[str, Any]] = []
    reports: List[Dict[str, Any]] = []

    try:
        if args.cmd == "daily":
            tag = f"recent_{args.recent_days}d"
            paths = debug_paths(tag)

            uids, esearch_url = esearch_recent("sra", args.query, args.recent_days, args.max_per_day, datetype="edat")
            summaries, esummary_url = esummary_sra(uids[:500])

            if args.debug:
                write_json(paths["initial"], {
                    "tag": tag, "query": args.query, "recent_days": args.recent_days,
                    "max_per_day": args.max_per_day, "esearch_url": esearch_url,
                    "esummary_url": esummary_url, "uids_count": len(uids), "uids_sample": uids[:50],
                })

            if uids:
                added_srr, report = ingest_uids_to_srr(
                    tag=tag, uids=uids, summaries=summaries,
                    biosample_cache=biosample_cache, bp_cache=bp_cache, bp_uid_cache=bp_uid_cache,
                    seen_sra=seen_sra, seen_srr=seen_srr,
                    fetch_biosample=args.fetch_biosample, fetch_bioproject=args.fetch_bioproject,
                    debug=args.debug, runinfo_max_rows=args.runinfo_max_rows,
                )
                if added_srr:
                    year = dt.date.today().year
                    append_jsonl(f"{DATA_DIR}/srr_catalog_{year}.jsonl", added_srr)

                    for r in added_srr[:800]:
                        bp = r.get("bioproject", {}) if isinstance(r.get("bioproject", {}), dict) else {}
                        latest_added.append({
                            "tag": tag,
                            "srr": r.get("srr", ""),
                            "sra_uid": r.get("sra_uid", ""),
                            "title": r.get("title", ""),
                            "assay_class": (r.get("assay") or {}).get("assay_class", ""),
                            "country": (r.get("geo") or {}).get("country", ""),
                            "city": (r.get("geo") or {}).get("city", ""),
                            "bioproject": (r.get("runinfo_row") or {}).get("BioProject", ""),
                            "bioproject_title": bp.get("title", ""),
                            "url": (r.get("ncbi") or {}).get("srr_url", ""),
                        })
                reports.append(report)
            else:
                reports.append({
                    "tag": tag,
                    "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                    "counters": {"uids_input": 0, "uids_new": 0, "srr_emitted": 0},
                    "urls": {"esearch": esearch_url, "esummary": esummary_url},
                })

        elif args.cmd == "backfill-year":
            year = args.year
            start = dt.date(year, 1, 1)
            for n in range(366):
                day = start + dt.timedelta(n)
                if day.year != year:
                    break
                ds = day.isoformat()

                uids, _ = esearch_day("sra", args.query, ds, args.max_per_day, datetype="edat")
                summaries, _ = esummary_sra(uids[:500])

                added_srr, report = ingest_uids_to_srr(
                    tag=ds, uids=uids, summaries=summaries,
                    biosample_cache=biosample_cache, bp_cache=bp_cache, bp_uid_cache=bp_uid_cache,
                    seen_sra=seen_sra, seen_srr=seen_srr,
                    fetch_biosample=args.fetch_biosample, fetch_bioproject=args.fetch_bioproject,
                    debug=args.debug, runinfo_max_rows=args.runinfo_max_rows,
                )
                if added_srr:
                    append_jsonl(f"{DATA_DIR}/srr_catalog_{year}.jsonl", added_srr)
                reports.append(report)

        else:  # crawl
            print("CRAWLING SEARCHES")
            retstart = 0
            page = 0
            total_seen = 0
            max_total = int(args.max_total or 0)
            new_srr_total = 0

            while True:
                ids, count_total, _ = esearch_history(
                    "sra", args.query, retstart=retstart, retmax=args.page_size, sort=args.sort
                )
                if not ids:
                    break

                summaries, _ = esummary_sra(ids[:500])
                tag = f"crawl_{page:06d}"

                print("[INFO] calling def ingest_uids_to_srr, count total=", count_total, "total seen=", total_seen)
                added_srr, report = ingest_uids_to_srr(
                    tag=tag, 
                    uids=ids, 
                    summaries=summaries,
                    biosample_cache=biosample_cache, 
                    bp_cache=bp_cache, 
                    bp_uid_cache=bp_uid_cache,
                    seen_sra=seen_sra, 
                    seen_srr=seen_srr,
                    fetch_biosample=args.fetch_biosample, 
                    fetch_bioproject=args.fetch_bioproject,
                    debug=args.debug, runinfo_max_rows=args.runinfo_max_rows,
                )

                new_srr_count = report["counters"].get("srr_emitted", 0)
                new_srr_total += new_srr_count

                if added_srr:
                    this_year = dt.date.today().year
                    append_jsonl(f"{DATA_DIR}/srr_catalog_{this_year}.jsonl", added_srr)

                reports.append(report)
                total_seen += len(ids)

                if args.stop_after_new_srr and new_srr_total >= args.stop_after_new_srr:
                    break

                retstart += args.page_size
                page += 1
                if retstart >= count_total:
                    break
                if max_total and total_seen >= max_total:
                    break

        write_latest_srr_safe(latest_added[:5000])
        write_json(DOCS_LATEST_DEBUG, {
            "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            "reports": reports
        })

        rebuild_srr_exports_chunked()

    finally:
        write_json(BIOSAMPLE_CACHE, biosample_cache)
        write_json(BIOPROJECT_CACHE, bp_cache)
        write_json(BIOPROJECT_UID_CACHE, bp_uid_cache)
