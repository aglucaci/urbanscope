#!/usr/bin/env python3
"""
UrbanScope — Urban Metagenomics & Metatranscriptomics Radar

Append-only discovery + enrichment pipeline for UNIQUE BioProjects (PRJ*),
built to run entirely on GitHub Actions and publish a static, exploreable
database via GitHub Pages.

Author: Alexander G. Lucaci
"""

from __future__ import annotations
import argparse, csv, json, os, re, time, random, urllib.request, urllib.parse
import datetime as dt
import xml.etree.ElementTree as ET
from typing import Dict, List, Set, Tuple, Iterable, Any

# ------------------------
# Configuration
# ------------------------
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

DEFAULT_QUERY = (
    '(urban OR city OR cities OR metropolis OR megacity OR subway OR transit '
    'OR "built environment" OR wastewater OR sewage OR stormwater) '
    'AND (metagenomics OR metatranscriptomics OR "shotgun metagenomic" OR RNA-seq)'
)

DATA_DIR = "data"
DOCS_DIR = "docs"
DB_DIR = f"{DOCS_DIR}/db"
DB_YEAR_DIR = f"{DB_DIR}/year"
ARCHIVE_DIR = f"{DOCS_DIR}/archive"
ARCHIVE_CSV_DIR = f"{ARCHIVE_DIR}/csv"

SEEN_IDS = f"{DATA_DIR}/seen_ids.txt"
SEEN_PROJECTS = f"{DATA_DIR}/seen_bioprojects.txt"

CACHE_DIR = f"{DATA_DIR}/cache"
BIOPROJECT_CACHE = f"{CACHE_DIR}/bioproject.json"
BIOSAMPLE_CACHE = f"{CACHE_DIR}/biosample.json"

BIOPROJECT_RE = re.compile(r"\bPRJ(?:NA|EB|DB)\d+\b", re.I)

# ------------------------
# Utilities
# ------------------------
def ensure_dirs():
    for p in [DATA_DIR, DOCS_DIR, DB_DIR, DB_YEAR_DIR, ARCHIVE_DIR, ARCHIVE_CSV_DIR, CACHE_DIR]:
        os.makedirs(p, exist_ok=True)

def utc_now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")

def http_get(url: str, retries: int = 5) -> bytes:
    for i in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                return r.read()
        except Exception:
            time.sleep(0.8 * (2 ** i) + random.random())
    raise RuntimeError(f"HTTP failed: {url}")

def parse_xml(data: bytes) -> ET.Element:
    return ET.fromstring(data)

def load_set(path: str) -> Set[str]:
    if not os.path.exists(path): return set()
    return set(x.strip() for x in open(path) if x.strip())

def append_lines(path: str, vals: Iterable[str]):
    with open(path, "a") as f:
        for v in vals:
            f.write(v + "\n")

def append_jsonl(path: str, records: List[Dict]):
    with open(path, "a", encoding="utf-8") as f:
        for r in records:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")

def iter_jsonl(path: str):
    if not os.path.exists(path): return
    for ln in open(path, encoding="utf-8"):
        if ln.strip():
            yield json.loads(ln)

# ------------------------
# NCBI helpers
# ------------------------
def esearch(db, term, day, retmax):
    params = {
        "db": db, "term": term, "retmode": "xml",
        "mindate": day, "maxdate": day, "datetype": "pdat",
        "retmax": retmax
    }
    url = EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(params)
    root = parse_xml(http_get(url))
    return [x.text for x in root.findall(".//Id") if x.text]

def esummary_sra(uids: List[str]) -> Dict[str, Dict]:
    if not uids: return {}
    url = EUTILS + "esummary.fcgi?" + urllib.parse.urlencode(
        {"db": "sra", "id": ",".join(uids), "retmode": "xml"}
    )
    root = parse_xml(http_get(url))
    out = {}
    for d in root.findall(".//DocSum"):
        uid = d.findtext("Id")
        title = ""
        bioproject = ""
        for it in d.findall("Item"):
            if it.attrib.get("Name") == "Title":
                title = it.text or ""
            if "BioProject" in it.attrib.get("Name", ""):
                m = BIOPROJECT_RE.search(it.text or "")
                if m: bioproject = m.group(0).upper()
        out[uid] = {"uid": uid, "title": title.strip(), "bioproject": bioproject}
    return out

def efetch_runinfo(uid: str) -> Dict[str, Any]:
    url = EUTILS + "efetch.fcgi?" + urllib.parse.urlencode(
        {"db": "sra", "id": uid, "rettype": "runinfo", "retmode": "text"}
    )
    text = http_get(url).decode()
    rows = list(csv.DictReader(text.splitlines()))
    biosamples = sorted({r.get("BioSample") for r in rows if r.get("BioSample")})
    platforms = {}
    for r in rows:
        p = r.get("Platform")
        if p: platforms[p] = platforms.get(p, 0) + 1
    return {
        "rows": len(rows),
        "biosample_count": len(biosamples),
        "biosamples": biosamples[:200],
        "platform_counts": platforms,
        "preview": rows[:20],
    }

# ------------------------
# Core logic
# ------------------------
def run_day(day: str, query: str, retmax: int,
            seen_ids: Set[str], seen_projects: Set[str]) -> Tuple[List[Dict], Set[str], Set[str]]:

    added, new_ids, new_projects = [], set(), set()

    sra_ids = esearch("sra", query, day, retmax)
    summaries = esummary_sra(sra_ids)

    for uid, s in summaries.items():
        bp = s["bioproject"]
        if not bp or bp in seen_projects or bp in new_projects:
            continue
        if not s["title"]:
            continue

        runinfo = efetch_runinfo(uid)

        rec = {
            "bioproject": bp,
            "title": s["title"],
            "representative_sra": uid,
            "runinfo": runinfo,
            "provenance": {
                "source": "sra",
                "day": day,
                "ingested_utc": utc_now(),
            }
        }

        added.append(rec)
        new_projects.add(bp)
        new_ids.add(f"sra:{uid}")

    return added, new_ids, new_projects

# ------------------------
# Exports
# ------------------------
def rebuild_exports():
    years = sorted(
        int(f.split("_")[1].split(".")[0])
        for f in os.listdir(DATA_DIR)
        if f.startswith("catalog_")
    )

    all_records = []
    for y in years:
        all_records.extend(list(iter_jsonl(f"{DATA_DIR}/catalog_{y}.jsonl")))

    with open(f"{DB_DIR}/records.json", "w") as f:
        json.dump(all_records, f, indent=2)

    with open(f"{DB_DIR}/index.json", "w") as f:
        json.dump({
            "generated": utc_now(),
            "total": len(all_records),
            "years": years
        }, f, indent=2)

# ------------------------
# CLI
# ------------------------
def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("backfill-year")
    b.add_argument("--year", type=int, required=True)
    b.add_argument("--query", default=DEFAULT_QUERY)
    b.add_argument("--max-per-day", type=int, default=200)

    d = sub.add_parser("daily")
    d.add_argument("--days", type=int, default=1)
    d.add_argument("--query", default=DEFAULT_QUERY)
    d.add_argument("--max-per-day", type=int, default=200)

    args = ap.parse_args()
    ensure_dirs()

    seen_ids = load_set(SEEN_IDS)
    seen_projects = load_set(SEEN_PROJECTS)

    if args.cmd == "backfill-year":
        year = args.year
        for day in (dt.date(year, 1, 1) + dt.timedelta(n) for n in range(365)):
            ds = day.isoformat()
            added, ids, projs = run_day(ds, args.query, args.max_per_day, seen_ids, seen_projects)
            if added:
                append_jsonl(f"{DATA_DIR}/catalog_{year}.jsonl", added)
            append_lines(SEEN_IDS, ids)
            append_lines(SEEN_PROJECTS, projs)
            seen_ids |= ids
            seen_projects |= projs
            print(ds, len(added))

    else:
        end = dt.date.today()
        for i in range(args.days):
            ds = (end - dt.timedelta(i)).isoformat()
            added, ids, projs = run_day(ds, args.query, args.max_per_day, seen_ids, seen_projects)
            if added:
                append_jsonl(f"{DATA_DIR}/catalog_{end.year}.jsonl", added)
            append_lines(SEEN_IDS, ids)
            append_lines(SEEN_PROJECTS, projs)
            seen_ids |= ids
            seen_projects |= projs
            print(ds, len(added))

    rebuild_exports()

if __name__ == "__main__":
    main()
