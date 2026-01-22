#!/usr/bin/env python3
"""
UrbanScope — Urban Metagenomics & Metatranscriptomics Radar

Append-only discovery + enrichment pipeline for UNIQUE BioProjects (PRJ*),
built to run entirely on GitHub Actions and publish a static, exploreable
database via GitHub Pages.

Key guarantees:
- Discovers SRA records by day, then resolves BioProject reliably via:
  (1) elink(sra -> bioproject) + esummary(bioproject -> accession)
  (2) fallback: efetch(runinfo) BioProject column (PRJ* accession)
- Enriches each BioProject with stable metadata + linked PubMed IDs.
- Uses local JSON caches in data/cache/ to minimize NCBI calls.

Author: Alexander G. Lucaci
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import time
import random
import urllib.request
import urllib.parse
import datetime as dt
import xml.etree.ElementTree as ET
from typing import Dict, List, Set, Tuple, Iterable, Any, Optional

# ------------------------
# Configuration
# ------------------------
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
NCBI_API_KEY = os.getenv("NCBI_API_KEY", "").strip()
TOOL_NAME = os.getenv("NCBI_TOOL", "urbanscope-radar")
NCBI_EMAIL = os.getenv("NCBI_EMAIL", "")

DEFAULT_QUERY = (
    '(urban OR city OR cities OR metropolis OR megacity OR subway OR transit '
    'OR "built environment" OR wastewater OR sewage OR stormwater) '
    'AND (metagenomics OR metatranscriptomics OR "shotgun metagenomic" OR RNA-seq)'
)

DATA_DIR = "data"
DOCS_DIR = "docs"
DB_DIR = f"{DOCS_DIR}/db"
ARCHIVE_DIR = f"{DOCS_DIR}/archive"
ARCHIVE_CSV_DIR = f"{ARCHIVE_DIR}/csv"

SEEN_IDS = f"{DATA_DIR}/seen_ids.txt"
SEEN_PROJECTS = f"{DATA_DIR}/seen_bioprojects.txt"

CACHE_DIR = f"{DATA_DIR}/cache"
BIOPROJECT_CACHE = f"{CACHE_DIR}/bioproject.json"          # accession -> details dict
BIOPROJECT_UID_CACHE = f"{CACHE_DIR}/bioproject_uid.json"  # accession -> uid (numeric)

BIOPROJECT_RE = re.compile(r"\bPRJ(?:NA|EB|DB)\d+\b", re.I)

# ------------------------
# Utilities
# ------------------------
def ensure_dirs():
    for p in [DATA_DIR, DOCS_DIR, DB_DIR, ARCHIVE_DIR, ARCHIVE_CSV_DIR, CACHE_DIR]:
        os.makedirs(p, exist_ok=True)

def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")

def _sleep_backoff(i: int):
    time.sleep(0.6 * (2 ** i) + random.random() * 0.25)

def http_get(url: str, retries: int = 6) -> bytes:
    headers = {"User-Agent": f"{TOOL_NAME}/1.0 ({NCBI_EMAIL or 'no-email'})"}
    req = urllib.request.Request(url, headers=headers)
    for i in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=40) as r:
                return r.read()
        except Exception:
            _sleep_backoff(i)
    raise RuntimeError(f"HTTP failed: {url}")

def parse_xml(data: bytes) -> ET.Element:
    return ET.fromstring(data)

def load_set(path: str) -> Set[str]:
    if not os.path.exists(path):
        return set()
    return set(x.strip() for x in open(path, encoding="utf-8") if x.strip())

def append_lines(path: str, vals: Iterable[str]):
    vals = list(vals)
    if not vals:
        return
    with open(path, "a", encoding="utf-8") as f:
        for v in vals:
            f.write(v + "\n")

def append_jsonl(path: str, records: List[Dict]):
    if not records:
        return
    with open(path, "a", encoding="utf-8") as f:
        for r in records:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")

def iter_jsonl(path: str):
    if not os.path.exists(path):
        return
    with open(path, encoding="utf-8") as fh:
        for ln in fh:
            if ln.strip():
                yield json.loads(ln)

def read_json(path: str, default):
    if not os.path.exists(path):
        return default
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return default

def write_json(path: str, obj):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)
    os.replace(tmp, path)

def chunked(xs: List[str], n: int) -> Iterable[List[str]]:
    for i in range(0, len(xs), n):
        yield xs[i:i+n]

# ------------------------
# NCBI helpers (E-utilities)
# ------------------------
def _eutils_params(extra: Dict[str, str]) -> Dict[str, str]:
    p = dict(extra)
    if NCBI_API_KEY:
        p["api_key"] = NCBI_API_KEY
    if TOOL_NAME:
        p["tool"] = TOOL_NAME
    if NCBI_EMAIL:
        p["email"] = NCBI_EMAIL
    return p

def esearch(db: str, term: str, day: str, retmax: int, datetype: str = "edat") -> List[str]:
    """
    SRA is often more consistent with edat than pdat.
    """
    params = _eutils_params({
        "db": db,
        "term": term,
        "retmode": "xml",
        "mindate": day,
        "maxdate": day,
        "datetype": datetype,
        "retmax": str(retmax),
    })
    url = EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(params)
    root = parse_xml(http_get(url))
    return [x.text for x in root.findall(".//Id") if x.text]

def esearch_any(db: str, term: str, retmax: int = 20) -> List[str]:
    params = _eutils_params({
        "db": db,
        "term": term,
        "retmode": "xml",
        "retmax": str(retmax),
    })
    url = EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(params)
    root = parse_xml(http_get(url))
    return [x.text for x in root.findall(".//Id") if x.text]

def esummary(db: str, ids: List[str]) -> ET.Element:
    if not ids:
        return ET.Element("EMPTY")
    params = _eutils_params({"db": db, "id": ",".join(ids), "retmode": "xml"})
    url = EUTILS + "esummary.fcgi?" + urllib.parse.urlencode(params)
    return parse_xml(http_get(url))

def elink(dbfrom: str, db: str, ids: List[str], linkname: Optional[str] = None) -> ET.Element:
    if not ids:
        return ET.Element("EMPTY")
    params: Dict[str, str] = {"dbfrom": dbfrom, "db": db, "id": ",".join(ids), "retmode": "xml"}
    if linkname:
        params["linkname"] = linkname
    params = _eutils_params(params)
    url = EUTILS + "elink.fcgi?" + urllib.parse.urlencode(params)
    return parse_xml(http_get(url))

def esummary_sra(uids: List[str]) -> Dict[str, Dict]:
    """
    Return per-SRA UID: title (BioProject accession in esummary is unreliable; we still attempt).
    """
    if not uids:
        return {}
    root = esummary("sra", uids)
    out: Dict[str, Dict] = {}
    for d in root.findall(".//DocSum"):
        uid = (d.findtext("Id") or "").strip()
        title = ""
        bioproject = ""  # best-effort only
        for it in d.findall("Item"):
            name = it.attrib.get("Name", "")
            if name == "Title":
                title = (it.text or "").strip()
            if "BioProject" in name:
                m = BIOPROJECT_RE.search(it.text or "")
                if m:
                    bioproject = m.group(0).upper()
        if uid:
            out[uid] = {"uid": uid, "title": title, "bioproject": bioproject}
    return out

def efetch_runinfo(uid: str) -> Dict[str, Any]:
    url = EUTILS + "efetch.fcgi?" + urllib.parse.urlencode(
        _eutils_params({"db": "sra", "id": uid, "rettype": "runinfo", "retmode": "text"})
    )
    text = http_get(url).decode(errors="replace")
    rows = list(csv.DictReader(text.splitlines()))

    biosamples = sorted({r.get("BioSample") for r in rows if r.get("BioSample")})
    platforms: Dict[str, int] = {}
    lib_strategies: Dict[str, int] = {}

    for r in rows:
        p = r.get("Platform")
        if p:
            platforms[p] = platforms.get(p, 0) + 1
        ls = r.get("LibraryStrategy")
        if ls:
            lib_strategies[ls] = lib_strategies.get(ls, 0) + 1

    keep_cols = [
        "Run", "BioProject", "BioSample", "SampleName",
        "LibraryStrategy", "LibrarySource", "LibrarySelection",
        "Platform", "Model", "spots", "bases", "avgLength",
        "ReleaseDate",
    ]
    preview = [{k: (r.get(k, "") or "") for k in keep_cols} for r in rows[:25]]

    return {
        "rows": len(rows),
        "biosample_count": len(biosamples),
        "biosamples": biosamples[:200],
        "platform_counts": platforms,
        "library_strategy_counts": lib_strategies,
        "preview": preview,
    }

# ------------------------
# Robust BioProject resolution from SRA UIDs
# ------------------------
def sra_uids_to_bioproject_uids(sra_uids: List[str]) -> Dict[str, str]:
    """
    Map SRA UID -> BioProject UID using elink. Returns {sra_uid: bioproject_uid}.
    """
    if not sra_uids:
        return {}

    out: Dict[str, str] = {}
    # elink can take many IDs; still chunk to be safe.
    for chunk in chunked(sra_uids, 200):
        root = elink("sra", "bioproject", chunk)
        for linkset in root.findall(".//LinkSet"):
            sra_id = (linkset.findtext("./IdList/Id") or "").strip()
            bp_id = (linkset.findtext(".//LinkSetDb/Link/Id") or "").strip()
            if sra_id and bp_id:
                out[sra_id] = bp_id
    return out

def bioproject_uid_to_accession(uid: str) -> str:
    """
    BioProject UID -> accession (PRJNA/PRJEB/PRJDB) via esummary.
    """
    uid = (uid or "").strip()
    if not uid:
        return ""
    root = esummary("bioproject", [uid])
    doc = root.find(".//DocSum")
    if doc is None:
        return ""
    for it in doc.findall("Item"):
        name = it.attrib.get("Name", "")
        if name in ("Project_Acc", "Accession"):
            acc = (it.text or "").strip()
            if acc:
                return acc.upper()
    return ""

def bioproject_from_runinfo(sra_uid: str) -> str:
    """
    Fallback: parse PRJ accession from SRA runinfo CSV's BioProject column.
    """
    url = EUTILS + "efetch.fcgi?" + urllib.parse.urlencode(
        _eutils_params({"db": "sra", "id": sra_uid, "rettype": "runinfo", "retmode": "text"})
    )
    text = http_get(url).decode(errors="replace")
    rows = list(csv.DictReader(text.splitlines()))
    for r in rows:
        m = BIOPROJECT_RE.search((r.get("BioProject") or ""))
        if m:
            return m.group(0).upper()
    return ""

# ------------------------
# BioProject enrichment (cached)
# ------------------------
def _load_bioproject_cache() -> Dict[str, Any]:
    return read_json(BIOPROJECT_CACHE, {})

def _load_bioproject_uid_cache() -> Dict[str, str]:
    return read_json(BIOPROJECT_UID_CACHE, {})

def _save_bioproject_cache(cache: Dict[str, Any]):
    write_json(BIOPROJECT_CACHE, cache)

def _save_bioproject_uid_cache(cache: Dict[str, str]):
    write_json(BIOPROJECT_UID_CACHE, cache)

def bioproject_accession_to_uid(accession: str, uid_cache: Dict[str, str]) -> Optional[str]:
    """
    Resolve PRJ accession -> BioProject numeric UID via esearch, cached.
    """
    accession = accession.upper()
    if accession in uid_cache:
        return uid_cache[accession] or None

    ids = esearch_any("bioproject", f"{accession}[Accession]", retmax=5)
    uid = ids[0] if ids else None
    uid_cache[accession] = uid or ""
    return uid

def parse_bioproject_esummary_docsum(docsum: ET.Element) -> Dict[str, Any]:
    uid = (docsum.findtext("Id") or "").strip()

    items: Dict[str, str] = {}
    for it in docsum.findall("Item"):
        name = it.attrib.get("Name", "")
        val = it.text or ""
        if list(it):
            subvals = [x.text for x in it.findall(".//Item") if x.text]
            if subvals:
                val = "; ".join(subvals)
        items[name] = val

    acc = (items.get("Project_Acc") or items.get("Accession") or "").strip()
    title = (items.get("Project_Title") or items.get("Title") or "").strip()
    desc = (items.get("Project_Description") or items.get("Description") or "").strip()
    organism = (items.get("Organism_Name") or items.get("Organism") or "").strip()
    data_type = (items.get("Project_Data_Type") or items.get("DataType") or "").strip()
    submission_date = (items.get("Submission_Date") or items.get("CreateDate") or "").strip()
    last_update = (items.get("Last_Update") or items.get("UpdateDate") or "").strip()
    center = (items.get("Center_Name") or items.get("Center") or items.get("Submitter") or "").strip()

    return {
        "uid": uid,
        "accession": acc,
        "title": title,
        "description": desc,
        "organism": organism,
        "data_type": data_type,
        "submission_date": submission_date,
        "last_update": last_update,
        "center_name": center,
        "ncbi": {
            "bioproject_uid": uid,
            "bioproject_url": f"https://www.ncbi.nlm.nih.gov/bioproject/{uid}" if uid else "",
        },
        "raw_esummary_keys": sorted(items.keys()),
    }

def fetch_bioproject_details(accession: str,
                             bp_cache: Dict[str, Any],
                             uid_cache: Dict[str, str],
                             pubmed: bool = True) -> Dict[str, Any]:
    accession = accession.upper()
    if accession in bp_cache:
        return bp_cache[accession]

    uid = bioproject_accession_to_uid(accession, uid_cache)
    details: Dict[str, Any] = {
        "uid": uid or "",
        "accession": accession,
        "title": "",
        "description": "",
        "organism": "",
        "data_type": "",
        "submission_date": "",
        "last_update": "",
        "center_name": "",
        "ncbi": {
            "bioproject_uid": uid or "",
            "bioproject_url": f"https://www.ncbi.nlm.nih.gov/bioproject/{uid}" if uid else "",
        },
        "linked_pubmed_ids": [],
    }

    if not uid:
        bp_cache[accession] = details
        return details

    root = esummary("bioproject", [uid])
    docsum = root.find(".//DocSum")
    if docsum is not None:
        parsed = parse_bioproject_esummary_docsum(docsum)
        parsed["accession"] = parsed.get("accession") or accession
        details.update(parsed)

    if pubmed:
        try:
            link_root = elink("bioproject", "pubmed", [uid])
            pmids = sorted({x.text for x in link_root.findall(".//LinkSetDb/Link/Id") if x.text})
            details["linked_pubmed_ids"] = pmids[:200]
        except Exception:
            details["linked_pubmed_ids"] = []

    bp_cache[accession] = details
    return details

# ------------------------
# Core logic
# ------------------------
def run_day(day: str,
            query: str,
            retmax: int,
            seen_ids: Set[str],
            seen_projects: Set[str],
            bp_cache: Dict[str, Any],
            uid_cache: Dict[str, str]) -> Tuple[List[Dict], Set[str], Set[str]]:

    added: List[Dict] = []
    new_ids: Set[str] = set()
    new_projects: Set[str] = set()

    sra_ids = esearch("sra", query, day, retmax, datetype="edat")
    if not sra_ids:
        return added, new_ids, new_projects

    summaries = esummary_sra(sra_ids)
    if not summaries:
        return added, new_ids, new_projects

    # Robust linkage: SRA UID -> BioProject UID
    sra_to_bpuid = sra_uids_to_bioproject_uids(list(summaries.keys()))

    # Convert unique bp uids -> PRJ accession
    bpuid_to_acc: Dict[str, str] = {}
    for bpuid in set(sra_to_bpuid.values()):
        acc = bioproject_uid_to_accession(bpuid)
        if acc:
            bpuid_to_acc[bpuid] = acc

    for uid, s in summaries.items():
        # Try (best-effort): esummary-supplied PRJ accession
        bp = (s.get("bioproject") or "").upper()

        # Try: elink-derived BioProject UID -> accession
        if not bp:
            bpuid = sra_to_bpuid.get(uid, "")
            bp = bpuid_to_acc.get(bpuid, "").upper()

        # Fallback: parse PRJ from runinfo CSV
        if not bp:
            bp = bioproject_from_runinfo(uid)

        if not bp:
            continue
        if bp in seen_projects or bp in new_projects:
            continue
        if not (s.get("title") or "").strip():
            continue

        runinfo = efetch_runinfo(uid)
        bp_details = fetch_bioproject_details(bp, bp_cache, uid_cache, pubmed=True)

        rec = {
            "bioproject": {
                "accession": bp,
                "uid": bp_details.get("uid", ""),
                "title": bp_details.get("title", "") or s["title"],
                "description": bp_details.get("description", ""),
                "organism": bp_details.get("organism", ""),
                "data_type": bp_details.get("data_type", ""),
                "submission_date": bp_details.get("submission_date", ""),
                "last_update": bp_details.get("last_update", ""),
                "center_name": bp_details.get("center_name", ""),
                "linked_pubmed_ids": bp_details.get("linked_pubmed_ids", []),
                "url": bp_details.get("ncbi", {}).get("bioproject_url", ""),
            },
            "representative_sra": {
                "uid": uid,
                "title": s["title"],
                "url": f"https://www.ncbi.nlm.nih.gov/sra/?term={uid}",
            },
            "run_summary": runinfo,
            "provenance": {
                "source": "sra",
                "day": day,
                "ingested_utc": utc_now(),
                "query": query,
            },
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
        if f.startswith("catalog_") and f.endswith(".jsonl")
    )

    all_records: List[Dict[str, Any]] = []
    for y in years:
        all_records.extend(list(iter_jsonl(f"{DATA_DIR}/catalog_{y}.jsonl")))

    with open(f"{DB_DIR}/records.json", "w", encoding="utf-8") as f:
        json.dump(all_records, f, ensure_ascii=False, indent=2)

    bp_map: Dict[str, Any] = {}
    for r in all_records:
        bp = r.get("bioproject")
        if isinstance(bp, dict):
            acc = (bp.get("accession") or "").upper()
            if acc and acc not in bp_map:
                bp_map[acc] = bp

    with open(f"{DB_DIR}/bioprojects.json", "w", encoding="utf-8") as f:
        json.dump(bp_map, f, ensure_ascii=False, indent=2)

    with open(f"{DB_DIR}/index.json", "w", encoding="utf-8") as f:
        json.dump(
            {
                "generated": utc_now(),
                "total_records": len(all_records),
                "total_bioprojects": len(bp_map),
                "years": years,
            },
            f,
            ensure_ascii=False,
            indent=2,
        )

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

    bp_cache = _load_bioproject_cache()
    uid_cache = _load_bioproject_uid_cache()

    try:
        if args.cmd == "backfill-year":
            year = args.year
            # Simple 365-day iterator (keeps your original behavior)
            for day in (dt.date(year, 1, 1) + dt.timedelta(n) for n in range(365)):
                ds = day.isoformat()
                added, ids, projs = run_day(ds, args.query, args.max_per_day, seen_ids, seen_projects, bp_cache, uid_cache)
                if added:
                    append_jsonl(f"{DATA_DIR}/catalog_{year}.jsonl", added)
                append_lines(SEEN_IDS, ids)
                append_lines(SEEN_PROJECTS, projs)
                seen_ids |= ids
                seen_projects |= projs
                print(ds, f"added={len(added)}")

        else:
            end = dt.date.today()
            for i in range(args.days):
                ds = (end - dt.timedelta(i)).isoformat()
                added, ids, projs = run_day(ds, args.query, args.max_per_day, seen_ids, seen_projects, bp_cache, uid_cache)
                if added:
                    append_jsonl(f"{DATA_DIR}/catalog_{end.year}.jsonl", added)
                append_lines(SEEN_IDS, ids)
                append_lines(SEEN_PROJECTS, projs)
                seen_ids |= ids
                seen_projects |= projs
                print(ds, f"added={len(added)}")

        rebuild_exports()

    finally:
        _save_bioproject_cache(bp_cache)
        _save_bioproject_uid_cache(uid_cache)

if __name__ == "__main__":
    main()
