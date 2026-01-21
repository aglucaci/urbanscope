#!/usr/bin/env python3
"""
UrbanScope Radar — historical backfill + daily incremental updates

Finds urban metagenomics/metatranscriptomics-related datasets by:
1) Searching SRA directly (db=sra) for the query in a date window
2) Searching PubMed for the query in a date window, then linking PubMed -> SRA via elink

It maintains:
- data/catalog.jsonl : append-only cumulative records
- data/seen_ids.txt  : dedupe index (SRA UIDs + PubMed IDs)

It writes:
- docs/latest.json : most recent items added in this run

Usage examples:
  # One-time historical backfill (2018-01-01 to 2025-12-31)
  python scripts/urbanscope_radar.py backfill --start 2018-01-01 --end 2025-12-31 --max-per-day 200

  # Daily incremental update (yesterday->today UTC)
  python scripts/urbanscope_radar.py daily --days 1 --max-per-day 200
"""

from __future__ import annotations
import argparse
import datetime as dt
import json
import os
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from typing import Dict, List, Any, Set, Tuple

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
DEFAULT_QUERY = (
    '(urban OR city OR cities OR subway OR transit OR "built environment" OR wastewater OR sewage) '
    'AND (metagenomics OR metatranscriptomics OR "shotgun metagenomic" OR "RNA-seq")'
)

DATA_DIR = "data"
DOCS_DIR = "docs"
SEEN_PATH = os.path.join(DATA_DIR, "seen_ids.txt")
CATALOG_PATH = os.path.join(DATA_DIR, "catalog.jsonl")
LATEST_PATH = os.path.join(DOCS_DIR, "latest.json")


def http_get(url: str, timeout: int = 30) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "urbanscope/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def esearch(db: str, term: str, mindate: str, maxdate: str, retmax: int, sort: str = "Date") -> List[str]:
    params = {
        "db": db,
        "term": term,
        "retmode": "xml",
        "retmax": str(retmax),
        "sort": sort,
        "mindate": mindate,
        "maxdate": maxdate,
        "datetype": "pdat",
    }
    url = EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(params)
    root = ET.fromstring(http_get(url))
    return [n.text for n in root.findall(".//IdList/Id") if n.text]


def elink_pubmed_to_sra(pmids: List[str]) -> Dict[str, List[str]]:
    if not pmids:
        return {}
    params = {
        "dbfrom": "pubmed",
        "db": "sra",
        "id": ",".join(pmids),
        "retmode": "xml",
    }
    url = EUTILS + "elink.fcgi?" + urllib.parse.urlencode(params)
    root = ET.fromstring(http_get(url))

    out: Dict[str, List[str]] = {}
    for linkset in root.findall("LinkSet"):
        pmid = (linkset.findtext("IdList/Id") or "").strip()
        sra_ids = [x.text.strip() for x in linkset.findall(".//Link/Id") if x.text]
        if pmid and sra_ids:
            out[pmid] = sra_ids
    return out


def esummary_sra(uids: List[str]) -> List[Dict[str, Any]]:
    """Pulls basic SRA docsum fields. We keep it lightweight."""
    if not uids:
        return []
    params = {"db": "sra", "id": ",".join(uids), "retmode": "xml"}
    url = EUTILS + "esummary.fcgi?" + urllib.parse.urlencode(params)
    root = ET.fromstring(http_get(url))
    items: List[Dict[str, Any]] = []

    for docsum in root.findall(".//DocSum"):
        uid = (docsum.findtext("Id") or "").strip()
        rec: Dict[str, Any] = {"sra_uid": uid}
        for it in docsum.findall("Item"):
            name = it.attrib.get("Name", "")
            text = (it.text or "").strip()
            # Keep a few helpful fields; ExpXml/Runs can be large so store only if present.
            if name in {"Title", "Runs"}:
                rec[name.lower()] = text
        rec["sra_link"] = f"https://www.ncbi.nlm.nih.gov/sra/?term={uid}" if uid else ""
        items.append(rec)
    return items


def load_seen(path: str) -> Set[str]:
    if not os.path.exists(path):
        return set()
    with open(path, "r", encoding="utf-8") as f:
        return set(line.strip() for line in f if line.strip())


def append_seen(path: str, ids: List[str]) -> None:
    if not ids:
        return
    with open(path, "a", encoding="utf-8") as f:
        for i in ids:
            f.write(i + "\n")


def append_catalog(path: str, records: List[Dict[str, Any]]) -> None:
    if not records:
        return
    with open(path, "a", encoding="utf-8") as f:
        for r in records:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")


def daterange(start: dt.date, end: dt.date) -> List[dt.date]:
    # inclusive start, inclusive end
    days = []
    cur = start
    while cur <= end:
        days.append(cur)
        cur += dt.timedelta(days=1)
    return days


def ensure_dirs() -> None:
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(DOCS_DIR, exist_ok=True)


def run_window(query: str, day_start: dt.date, day_end: dt.date, retmax: int, seen: Set[str]) -> Tuple[List[Dict[str, Any]], List[str]]:
    """
    For a given date window, collect:
    - SRA hits directly
    - PubMed hits linked to SRA
    Return new catalog records + newly seen ids to add to seen file.
    """
    mindate = day_start.strftime("%Y/%m/%d")
    maxdate = day_end.strftime("%Y/%m/%d")

    new_records: List[Dict[str, Any]] = []
    new_seen: List[str] = []

    # A) Search SRA directly
    try:
        sra_uids = esearch("sra", query, mindate, maxdate, retmax=retmax, sort="Date")
    except Exception:
        sra_uids = []
    time.sleep(0.34)

    sra_uids_new = [u for u in sra_uids if f"sra:{u}" not in seen]
    if sra_uids_new:
        summaries = esummary_sra(sra_uids_new[:200])  # keep API payload reasonable
        for s in summaries:
            uid = s.get("sra_uid", "")
            if not uid:
                continue
            rec = {
                "source": "sra_search",
                "date_window": {"mindate": mindate, "maxdate": maxdate},
                "query": query,
                "sra_uid": uid,
                "sra_link": s.get("sra_link", ""),
                "title": s.get("title", ""),
                "runs": s.get("runs", ""),
                "ingested_utc": dt.datetime.utcnow().isoformat(timespec="seconds"),
            }
            new_records.append(rec)
            new_seen.append(f"sra:{uid}")
    time.sleep(0.34)

    # B) Search PubMed then link to SRA
    try:
        pmids = esearch("pubmed", query, mindate, maxdate, retmax=retmax, sort="pub+date")
    except Exception:
        pmids = []
    time.sleep(0.34)

    pmids_new = [p for p in pmids if f"pmid:{p}" not in seen]
    if pmids_new:
        pmid_to_sra = elink_pubmed_to_sra(pmids_new[:200])
        for pmid, uids in pmid_to_sra.items():
            # Record the linkage; datasets are SRA uids
            rec = {
                "source": "pubmed_elink_sra",
                "date_window": {"mindate": mindate, "maxdate": maxdate},
                "query": query,
                "pmid": pmid,
                "paper": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                "sra_uids": uids,
                "sra_link": f"https://www.ncbi.nlm.nih.gov/sra/?term={' OR '.join(uids)}" if uids else "",
                "ingested_utc": dt.datetime.utcnow().isoformat(timespec="seconds"),
            }
            new_records.append(rec)
            new_seen.append(f"pmid:{pmid}")
            for uid in uids:
                new_seen.append(f"sra:{uid}")
    time.sleep(0.34)

    # Deduplicate within-run
    dedup_seen = []
    seen_local = set()
    for x in new_seen:
        if x not in seen and x not in seen_local:
            dedup_seen.append(x)
            seen_local.add(x)

    # Also ensure we don't emit duplicate records for already-seen datasets
    filtered_records = []
    for r in new_records:
        # Basic filter: if record is a direct SRA record, ensure sra uid not already seen
        if r.get("source") == "sra_search":
            uid = r.get("sra_uid", "")
            if uid and f"sra:{uid}" in seen:
                continue
        filtered_records.append(r)

    return filtered_records, dedup_seen


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_back = sub.add_parser("backfill", help="Historical backfill (date range), appended into catalog")
    p_back.add_argument("--start", required=True, help="YYYY-MM-DD")
    p_back.add_argument("--end", required=True, help="YYYY-MM-DD")
    p_back.add_argument("--query", default=DEFAULT_QUERY)
    p_back.add_argument("--max-per-day", type=int, default=200)

    p_daily = sub.add_parser("daily", help="Daily incremental update (last N days)")
    p_daily.add_argument("--days", type=int, default=1)
    p_daily.add_argument("--query", default=DEFAULT_QUERY)
    p_daily.add_argument("--max-per-day", type=int, default=200)

    args = ap.parse_args()
    ensure_dirs()

    seen = load_seen(SEEN_PATH)
    added_records_all: List[Dict[str, Any]] = []
    added_ids_all: List[str] = []

    if args.cmd == "backfill":
        start = dt.date.fromisoformat(args.start)
        end = dt.date.fromisoformat(args.end)
        days = daterange(start, end)
        # Walk day-by-day so the search window is stable and easy to resume.
        for d in days:
            recs, new_seen = run_window(args.query, d, d, args.max_per_day, seen)
            if recs:
                append_catalog(CATALOG_PATH, recs)
                added_records_all.extend(recs)
            if new_seen:
                append_seen(SEEN_PATH, new_seen)
                for s in new_seen:
                    seen.add(s)
                added_ids_all.extend(new_seen)
            # Progress print
            print(f"[{d.isoformat()}] +{len(recs)} records, +{len(new_seen)} ids")
    else:
        # daily: last N days up to today UTC inclusive
        end = dt.datetime.utcnow().date()
        start = end - dt.timedelta(days=max(args.days, 1))
        for d in daterange(start, end):
            recs, new_seen = run_window(args.query, d, d, args.max_per_day, seen)
            if recs:
                append_catalog(CATALOG_PATH, recs)
                added_records_all.extend(recs)
            if new_seen:
                append_seen(SEEN_PATH, new_seen)
                for s in new_seen:
                    seen.add(s)
                added_ids_all.extend(new_seen)
            print(f"[{d.isoformat()}] +{len(recs)} records, +{len(new_seen)} ids")

    # Write latest.json (what was added this run)
    latest = {
        "generated_at_utc": dt.datetime.utcnow().isoformat(timespec="seconds"),
        "query": args.query,
        "added_records": len(added_records_all),
        "added_ids": len(added_ids_all),
        "records": added_records_all[-200:],  # keep latest.json small
    }
    with open(LATEST_PATH, "w", encoding="utf-8") as f:
        json.dump(latest, f, indent=2, ensure_ascii=False)

    print(f"Wrote {LATEST_PATH} (added_records={len(added_records_all)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
