from __future__ import annotations
import argparse, datetime as dt, json
from typing import Any, Dict, List

from .config import (
    DEFAULT_QUERY, QUERY_PROFILES, DEFAULT_QUERY_PROFILE_NAMES,
    BIOSAMPLE_CACHE, BIOPROJECT_CACHE, BIOPROJECT_UID_CACHE, AI_CURATION_CACHE, OPENAI_MODEL, OPENAI_API_KEY,
    DOCS_LATEST_DEBUG, DATA_DIR, DB_DIR
)
from .utils import ensure_dirs, load_set, read_json, write_json, append_jsonl, iter_jsonl_glob
from .ncbi import esearch_recent, esearch_day, esearch_history, esummary_sra
from .ingest import ingest_uids_to_srr, debug_paths
from .exports import rebuild_srr_exports_chunked, write_latest_srr_safe
from .ai_curation import curate_records

def print_report_summary(report: Dict[str, Any]):
    tag = report.get("tag", "run")
    counters = report.get("counters", {}) if isinstance(report.get("counters", {}), dict) else {}
    ai = report.get("ai_curation", {}) if isinstance(report.get("ai_curation", {}), dict) else {}
    parts = [f"tag={tag}"]

    if "records_considered" in report:
        parts.append(f"records_considered={report.get('records_considered', 0)}")
    if counters:
        parts.extend([
            f"uids_input={counters.get('uids_input', 0)}",
            f"uids_new={counters.get('uids_new', 0)}",
            f"srr_emitted={counters.get('srr_emitted', 0)}",
        ])
    if ai:
        parts.extend([
            f"ai_reviewed={ai.get('reviewed', 0)}",
            f"ai_skipped_cached={ai.get('skipped_cached', 0)}",
            f"ai_errors={ai.get('errors', 0)}",
        ])

    print("[SUMMARY] " + " ".join(parts))
    print("[SUMMARY_JSON] " + json.dumps(report, ensure_ascii=False))

def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    def add_ai_args(parser):
        parser.add_argument("--ai-curate", action="store_true")
        parser.add_argument("--ai-model", default=OPENAI_MODEL)
        parser.add_argument("--ai-max-records", type=int, default=0)

    d = sub.add_parser("daily")
    d.add_argument("--days", type=int, default=1)  # kept for compatibility
    d.add_argument("--query", default=DEFAULT_QUERY)
    d.add_argument("--query-profile", action="append", choices=sorted(QUERY_PROFILES.keys()))
    d.add_argument("--max-per-day", type=int, default=500)
    d.add_argument("--debug", action="store_true")
    d.add_argument("--fetch-biosample", action="store_true")
    d.add_argument("--fetch-bioproject", action="store_true")
    d.add_argument("--runinfo-max-rows", type=int, default=200000)
    d.add_argument("--recent-days", type=int, default=7)
    add_ai_args(d)

    b = sub.add_parser("backfill-year")
    b.add_argument("--year", type=int, required=True)
    b.add_argument("--query", default=DEFAULT_QUERY)
    b.add_argument("--query-profile", action="append", choices=sorted(QUERY_PROFILES.keys()))
    b.add_argument("--max-per-day", type=int, default=500)
    b.add_argument("--debug", action="store_true")
    b.add_argument("--fetch-biosample", action="store_true")
    b.add_argument("--fetch-bioproject", action="store_true")
    b.add_argument("--runinfo-max-rows", type=int, default=200000)
    add_ai_args(b)

    c = sub.add_parser("crawl")
    c.add_argument("--query", default=DEFAULT_QUERY)
    c.add_argument("--query-profile", action="append", choices=sorted(QUERY_PROFILES.keys()))
    c.add_argument("--page-size", type=int, default=500)
    c.add_argument("--max-total", type=int, default=0)
    c.add_argument("--debug", action="store_true")
    c.add_argument("--fetch-biosample", action="store_true")
    c.add_argument("--fetch-bioproject", action="store_true")
    c.add_argument("--runinfo-max-rows", type=int, default=200000)
    c.add_argument("--stop-after-new-srr", type=int, default=0)
    c.add_argument("--sort", default="date")
    add_ai_args(c)

    a = sub.add_parser("curate-ai")
    a.add_argument("--model", default=OPENAI_MODEL)
    a.add_argument("--max-records", type=int, default=0)
    a.add_argument("--overwrite", action="store_true")
    a.add_argument("--year", type=int, default=0)
    return ap

def resolve_query_specs(args) -> List[Dict[str, str]]:
    selected_profiles = getattr(args, "query_profile", None) or []
    query = getattr(args, "query", DEFAULT_QUERY)
    use_default_profiles = (query == DEFAULT_QUERY) and not selected_profiles

    if use_default_profiles:
      return [{"name": name, "query": QUERY_PROFILES[name]} for name in DEFAULT_QUERY_PROFILE_NAMES]

    if selected_profiles:
      return [{"name": name, "query": QUERY_PROFILES[name]} for name in selected_profiles]

    return [{"name": "custom", "query": query}]

def run():
    args = build_argparser().parse_args()
    ensure_dirs()

    if ((getattr(args, "ai_curate", False) or args.cmd == "curate-ai") and not OPENAI_API_KEY):
        raise RuntimeError("OPENAI_API_KEY is required for AI-assisted curation")

    from .config import SEEN_SRA_UIDS, SEEN_SRR_RUNS, DOCS_LATEST_SRR

    seen_sra = load_set(SEEN_SRA_UIDS)
    seen_srr = load_set(SEEN_SRR_RUNS)

    biosample_cache = read_json(BIOSAMPLE_CACHE, {})
    bp_cache = read_json(BIOPROJECT_CACHE, {})
    bp_uid_cache = read_json(BIOPROJECT_UID_CACHE, {})
    ai_cache = read_json(AI_CURATION_CACHE, {})

    latest_added: List[Dict[str, Any]] = []
    reports: List[Dict[str, Any]] = []
    query_specs = resolve_query_specs(args)

    try:
        if args.cmd == "daily":
            tag = f"recent_{args.recent_days}d"
            paths = debug_paths(tag)
            uids: List[str] = []
            query_reports: List[Dict[str, Any]] = []
            seen_uids = set()
            for spec in query_specs:
                found_uids, esearch_url = esearch_recent("sra", spec["query"], args.recent_days, args.max_per_day, datetype="edat")
                query_reports.append({
                    "profile": spec["name"],
                    "query": spec["query"],
                    "uids_count": len(found_uids),
                    "esearch_url": esearch_url,
                })
                for uid in found_uids:
                    if uid not in seen_uids:
                        seen_uids.add(uid)
                        uids.append(uid)
            summaries, esummary_url = esummary_sra(uids[:500])

            if args.debug:
                write_json(paths["initial"], {
                    "tag": tag, "queries": query_reports, "recent_days": args.recent_days,
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
                    if args.ai_curate:
                        ai_counts = curate_records(
                            added_srr,
                            ai_cache=ai_cache,
                            biosample_cache=biosample_cache,
                            bioproject_cache=bp_cache,
                            model=args.ai_model,
                            max_records=args.ai_max_records,
                        )
                        report.setdefault("ai_curation", ai_counts)
                    print_report_summary(report)
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
                report = {
                    "tag": tag,
                    "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                    "counters": {"uids_input": 0, "uids_new": 0, "srr_emitted": 0},
                    "queries": query_reports,
                    "urls": {"esummary": esummary_url},
                }
                print_report_summary(report)
                reports.append(report)

        elif args.cmd == "backfill-year":
            year = args.year
            start = dt.date(year, 1, 1)
            for n in range(366):
                day = start + dt.timedelta(n)
                if day.year != year:
                    break
                ds = day.isoformat()
                uids: List[str] = []
                seen_uids = set()
                for spec in query_specs:
                    found_uids, _ = esearch_day("sra", spec["query"], ds, args.max_per_day, datetype="edat")
                    for uid in found_uids:
                        if uid not in seen_uids:
                            seen_uids.add(uid)
                            uids.append(uid)
                summaries, _ = esummary_sra(uids[:500])

                added_srr, report = ingest_uids_to_srr(
                    tag=ds, uids=uids, summaries=summaries,
                    biosample_cache=biosample_cache, bp_cache=bp_cache, bp_uid_cache=bp_uid_cache,
                    seen_sra=seen_sra, seen_srr=seen_srr,
                    fetch_biosample=args.fetch_biosample, fetch_bioproject=args.fetch_bioproject,
                    debug=args.debug, runinfo_max_rows=args.runinfo_max_rows,
                )
                if added_srr:
                    if args.ai_curate:
                        ai_counts = curate_records(
                            added_srr,
                            ai_cache=ai_cache,
                            biosample_cache=biosample_cache,
                            bioproject_cache=bp_cache,
                            model=args.ai_model,
                            max_records=args.ai_max_records,
                        )
                        report.setdefault("ai_curation", ai_counts)
                    append_jsonl(f"{DATA_DIR}/srr_catalog_{year}.jsonl", added_srr)
                print_report_summary(report)
                reports.append(report)

        elif args.cmd == "crawl":
            print("CRAWLING SEARCHES")
            retstart = 0
            page = 0
            total_seen = 0
            max_total = int(args.max_total or 0)
            new_srr_total = 0

            while True:
                ids: List[str] = []
                count_total = 0
                seen_page_ids = set()
                for spec in query_specs:
                    found_ids, found_total, _ = esearch_history(
                        "sra", spec["query"], retstart=retstart, retmax=args.page_size, sort=args.sort
                    )
                    count_total = max(count_total, found_total)
                    for uid in found_ids:
                        if uid not in seen_page_ids:
                            seen_page_ids.add(uid)
                            ids.append(uid)
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
                    if args.ai_curate:
                        ai_counts = curate_records(
                            added_srr,
                            ai_cache=ai_cache,
                            biosample_cache=biosample_cache,
                            bioproject_cache=bp_cache,
                            model=args.ai_model,
                            max_records=args.ai_max_records,
                        )
                        report.setdefault("ai_curation", ai_counts)
                    this_year = dt.date.today().year
                    append_jsonl(f"{DATA_DIR}/srr_catalog_{this_year}.jsonl", added_srr)

                print_report_summary(report)
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

        else:  # curate-ai
            records: List[Dict[str, Any]] = []
            year_prefix = f"srr_catalog_{args.year}.jsonl" if args.year else None
            from .exports import _find_year_catalog_prefixes
            for base in _find_year_catalog_prefixes():
                if year_prefix and year_prefix not in base:
                    continue
                for rec in iter_jsonl_glob(base):
                    records.append(rec)
                    if args.max_records and len(records) >= args.max_records:
                        break
                if args.max_records and len(records) >= args.max_records:
                    break

            if not records:
                manifest = read_json(f"{DB_DIR}/srr_records_manifest.json", {})
                parts = manifest.get("parts", []) if isinstance(manifest, dict) else []
                for part in parts:
                    part_path = part.get("path", "") if isinstance(part, dict) else ""
                    if not part_path:
                        continue
                    chunk = read_json(part_path, [])
                    if not isinstance(chunk, list):
                        continue
                    for rec in chunk:
                        records.append(rec)
                        if args.max_records and len(records) >= args.max_records:
                            break
                    if args.max_records and len(records) >= args.max_records:
                        break

            ai_counts = curate_records(
                records,
                ai_cache=ai_cache,
                biosample_cache=biosample_cache,
                bioproject_cache=bp_cache,
                model=args.model,
                overwrite=args.overwrite,
                max_records=args.max_records,
            )
            report = {
                "tag": "curate_ai",
                "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                "ai_curation": ai_counts,
                "records_considered": len(records),
            }
            print_report_summary(report)
            reports.append(report)

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
        write_json(AI_CURATION_CACHE, ai_cache)
