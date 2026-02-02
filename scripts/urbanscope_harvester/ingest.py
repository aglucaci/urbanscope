from __future__ import annotations
import csv
from typing import Any, Dict, List, Set, Tuple

from .utils import inc, utc_now, append_lines, append_jsonl_one
from .config import SEEN_SRA_UIDS, SEEN_SRR_RUNS, BIOPROJECT_RE
from .ncbi import efetch_runinfo_text
from .biosample import get_biosample_details, infer_geo
from .assay import classify_assay
from .bioproject import get_bioproject_details

def parse_runinfo_rows(uid: str, max_rows: int = 200000) -> Tuple[List[Dict[str, str]], Dict[str, Any]]:
    text, url = efetch_runinfo_text(uid)
    rows = list(csv.DictReader(text.splitlines()))
    if max_rows and len(rows) > max_rows:
        rows = rows[:max_rows]
    cols = list(rows[0].keys()) if rows else []
    return rows, {"url": url, "columns": cols, "rows": len(rows)}

def debug_paths(tag: str) -> Dict[str, str]:
    from .config import DEBUG_DIR
    return {
        "report": f"{DEBUG_DIR}/report_{tag}.json",
        "decision": f"{DEBUG_DIR}/decision_log_{tag}.ndjson",
        "initial": f"{DEBUG_DIR}/initial_{tag}.json",
    }

def build_srr_records_for_sra_uid(
    sra_uid: str,
    sra_summary: Dict[str, Any],
    biosample_cache: Dict[str, Any],
    fetch_biosample: bool,
    bp_cache: Dict[str, Any],
    bp_uid_cache: Dict[str, str],
    fetch_bioproject: bool,
    debug: bool,
    decision_log_path: str,
    runinfo_max_rows: int,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    rows, runinfo_dbg = parse_runinfo_rows(sra_uid, max_rows=runinfo_max_rows)
    title = (sra_summary.get("title") or "").strip()
    #print("in build_srr_records_for_sra_uid")
    out: List[Dict[str, Any]] = []
    for r in rows:
        srr = (r.get("Run") or "").strip()
        if not srr:
            continue

        biosample_acc = (r.get("BioSample") or "").strip()
        biosample_details = {}
        if fetch_biosample and biosample_acc:
            biosample_details = get_biosample_details(biosample_acc, biosample_cache)

        geo = infer_geo(
            biosample_details,
            fallbacks=[
                title,
                r.get("SampleName", "") or "",
                r.get("Sample", "") or "",
                r.get("Study", "") or "",
                r.get("BioProject", "") or "",
            ],
        )

        assay = classify_assay(r, title, biosample_details)

        prj = (r.get("BioProject") or "").strip().upper()
        if not prj and (sra_summary.get("bioproject_guess") or ""):
            prj = (sra_summary.get("bioproject_guess") or "").strip().upper()

        bioproject_details = {}
        if fetch_bioproject and prj and BIOPROJECT_RE.match(prj):
            bioproject_details = get_bioproject_details(prj, bp_cache, bp_uid_cache)

        out.append({
            "srr": srr,
            "sra_uid": sra_uid,
            "title": title,
            "runinfo_row": r,
            "geo": geo,
            "assay": assay,
            "bioproject": bioproject_details,
            "ncbi": {
                "sra_uid_url": f"https://www.ncbi.nlm.nih.gov/sra/?term={sra_uid}",
                "srr_url": f"https://www.ncbi.nlm.nih.gov/sra/?term={srr}",
                "bioproject_url": bioproject_details.get("ncbi", {}).get("bioproject_url", "") if isinstance(bioproject_details, dict) else "",
            },
            "debug": {
                "runinfo_url": runinfo_dbg.get("url", ""),
                "runinfo_columns": runinfo_dbg.get("columns", []),
                "biosample_url": (biosample_details or {}).get("efetch_url", ""),
            },
            "provenance": {"ingested_utc": utc_now(), "source": "ncbi_eutils"},
        })

    if debug:
        append_jsonl_one(decision_log_path, {
            "uid": sra_uid,
            "decision": "flattened_to_srr",
            "runs_emitted": len(out),
            "runinfo_url": runinfo_dbg.get("url", ""),
            "runinfo_rows": runinfo_dbg.get("rows", 0),
        })

    return out, runinfo_dbg

def ingest_uids_to_srr(
    tag: str,
    uids: List[str],
    summaries: Dict[str, Dict[str, Any]],
    biosample_cache: Dict[str, Any],
    bp_cache: Dict[str, Any],
    bp_uid_cache: Dict[str, str],
    seen_sra: Set[str],
    seen_srr: Set[str],
    fetch_biosample: bool,
    fetch_bioproject: bool,
    debug: bool,
    runinfo_max_rows: int,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    from .utils import write_json
    counters: Dict[str, int] = {}
    paths = debug_paths(tag)
 
    #print("IN ingest_uids_to_srr")

    inc(counters, "uids_input", len(uids))
    new_uids = [u for u in uids if u not in seen_sra]
    inc(counters, "uids_new", len(new_uids))

    if debug:
        write_json(paths["initial"], {
            "tag": tag,
            "uids_input_count": len(uids),
            "uids_new_count": len(new_uids),
            "uids_new_sample": new_uids[:50],
        })

    added_srr: List[Dict[str, Any]] = []
    sra_mark_seen: List[str] = []
    srr_mark_seen: List[str] = []

    for uid in new_uids:
        ssum = summaries.get(uid, {"uid": uid, "title": "", "bioproject_guess": "", "items": {}})
        try:
            srr_rows, _ = build_srr_records_for_sra_uid(
                sra_uid=uid,
                sra_summary=ssum,
                biosample_cache=biosample_cache,
                fetch_biosample=fetch_biosample,
                bp_cache=bp_cache,
                bp_uid_cache=bp_uid_cache,
                fetch_bioproject=fetch_bioproject,
                debug=debug,
                decision_log_path=paths["decision"],
                runinfo_max_rows=runinfo_max_rows,
            )

            emitted = 0
            for r in srr_rows:
                srr = (r.get("srr") or "").strip()
                if not srr:
                    continue
                if srr in seen_srr:
                    inc(counters, "skip_seen_srr")
                    continue
                added_srr.append(r)
                seen_srr.add(srr)
                srr_mark_seen.append(srr)
                emitted += 1

            inc(counters, "srr_emitted", emitted)
            inc(counters, "sra_uids_processed_ok")

            seen_sra.add(uid)
            sra_mark_seen.append(uid)

        except Exception as e:
            inc(counters, "sra_uid_errors")
            if debug:
                append_jsonl_one(paths["decision"], {"uid": uid, "decision": "error", "error": str(e)})

    if sra_mark_seen:
        append_lines(SEEN_SRA_UIDS, sra_mark_seen)
    if srr_mark_seen:
        append_lines(SEEN_SRR_RUNS, srr_mark_seen)

    report = {"tag": tag, "generated_utc": utc_now(), "counters": counters, "debug_files": paths if debug else {}}
    if debug:
        write_json(paths["report"], report)

    return added_srr, report
