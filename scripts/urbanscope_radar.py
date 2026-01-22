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
DB_YEAR_DIR = f"{DB_DIR}/year"
ARCHIVE_DIR = f"{DOCS_DIR}/archive"
ARCHIVE_CSV_DIR = f"{ARCHIVE_DIR}/csv"

SEEN_IDS = f"{DATA_DIR}/seen_ids.txt"
SEEN_PROJECTS = f"{DATA_DIR}/seen_bioprojects.txt"

CACHE_DIR = f"{DATA_DIR}/cache"
BIOPROJECT_CACHE = f"{CACHE_DIR}/bioproject.json"   # accession -> details dict
BIOPROJECT_UID_CACHE = f"{CACHE_DIR}/bioproject_uid.json"  # accession -> uid (numeric)

BIOPROJECT_RE = re.compile(r"\bPRJ(?:NA|EB|DB)\d+\b", re.I)

# ------------------------
# Utilities
# ------------------------
def ensure_dirs():
    for p in [DATA_DIR, DOCS_DIR, DB_DIR, DB_YEAR_DIR, ARCHIVE_DIR, ARCHIVE_CSV_DIR, CACHE_DIR]:
        os.makedirs(p, exist_ok=True)

def utc_now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")

def _sleep_backoff(i: int):
    time.sleep(0.6 * (2 ** i) + random.random() * 0.25)

def http_get(url: str, retries: int = 6) -> bytes:
    # NCBI likes a User-Agent; include tool/email if available.
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
    # Optional but recommended by NCBI
    if NCBI_API_KEY:
        p["api_key"] = NCBI_API_KEY
    if TOOL_NAME:
        p["tool"] = TOOL_NAME
    if NCBI_EMAIL:
        p["email"] = NCBI_EMAIL
    return p

def esearch(db: str, term: str, day: str, retmax: int) -> List[str]:
    params = _eutils_params({
        "db": db,
        "term": term,
        "retmode": "xml",
        "mindate": day,
        "maxdate": day,
        "datetype": "pdat",
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
    params = {"dbfrom": dbfrom, "db": db, "id": ",".join(ids), "retmode": "xml"}
    if linkname:
        params["linkname"] = linkname
    params = _eutils_params({k: str(v) for k, v in params.items()})
    url = EUTILS + "elink.fcgi?" + urllib.parse.urlencode(params)
    return parse_xml(http_get(url))

def esummary_sra(uids: List[str]) -> Dict[str, Dict]:
    """
    Return per-SRA UID: title + first BioProject accession discovered in the summary.
    """
    if not uids:
        return {}
    root = esummary("sra", uids)
    out: Dict[str, Dict] = {}
    for d in root.findall(".//DocSum"):
        uid = d.findtext("Id") or ""
        title = ""
        bioproject = ""
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

    # Some useful “preview” columns (keeps the JSON smaller + stable for the UI)
    keep_cols = [
        "Run", "BioProject", "BioSample", "SampleName",
        "LibraryStrategy", "LibrarySource", "LibrarySelection",
        "Platform", "Model", "spots", "bases", "avgLength",
        "ReleaseDate",
    ]
    preview = []
    for r in rows[:25]:
        preview.append({k: r.get(k, "") for k in keep_cols})

    return {
        "rows": len(rows),
        "biosample_count": len(biosamples),
        "biosamples": biosamples[:200],
        "platform_counts": platforms,
        "library_strategy_counts": lib_strategies,
        "preview": preview,
    }

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

def bioproject_accession_to_uid(accession: str,
                                uid_cache: Dict[str, str]) -> Optional[str]:
    """
    BioProject esummary/efetch uses numeric UIDs. We resolve PRJNA... -> UID via esearch.
    Cached to avoid repeated lookups.
    """
    accession = accession.upper()
    if accession in uid_cache:
        return uid_cache[accession] or None

    # Use [Accession] query to resolve reliably
    term = f"{accession}[Accession]"
    ids = esearch_any("bioproject", term, retmax=5)
    uid = ids[0] if ids else None
    uid_cache[accession] = uid or ""
    return uid

def parse_bioproject_esummary_docsum(docsum: ET.Element) -> Dict[str, Any]:
    """
    Extract useful BioProject metadata from esummary DocSum.
    Fields vary a bit; we keep it resilient.
    """
    uid = docsum.findtext("Id") or ""
    items = {}
    for it in docsum.findall("Item"):
        name = it.attrib.get("Name", "")
        val = it.text or ""
        # Some values are nested lists; attempt to flatten a few cases
        if list(it):
            # join nested Items' text where applicable
            subvals = [x.text for x in it.findall(".//Item") if x.text]
            if subvals:
                val = "; ".join(subvals)
        items[name] = val

    # Commonly present in BioProject esummary:
    # Project_Acc, Project_Title, Project_Description, Organism_Name, Project_Data_Type,
    # Submission_Date, Last_Update, Center_Name, Accession (sometimes), etc.
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
    """
    Return a stable, UI-friendly BioProject details dict, cached by accession.
    Optionally add linked PubMed IDs (lightweight and useful in the UI).
    """
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

    # esummary gives the most consistent structured fields for BioProject
    root = esummary("bioproject", [uid])
    docsum = root.find(".//DocSum")
    if docsum is not None:
        parsed = parse_bioproject_esummary_docsum(docsum)
        # ensure accession field is set even if esummary omits it
        parsed["accession"] = parsed.get("accession") or accession
        details.update(parsed)

    # Optional: linked PubMed IDs (helps your database “evidence” layer)
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

    sra_ids = esearch("sra", query, day, retmax)
    summaries = esummary_sra(sra_ids)

    for uid, s in summaries.items():
        bp = (s.get("bioproject") or "").upper()
        if not bp:
            continue
        if bp in seen_projects or bp in new_projects:
            continue
        if not (s.get("title") or "").strip():
            continue

        # Enrich: run-level overview + BioProject-level metadata
        runinfo = efetch_runinfo(uid)
        bp_details = fetch_bioproject_details(bp, bp_cache, uid_cache, pubmed=True)

        rec = {
            # Primary key for your database
            "bioproject": {
                "accession": bp,
                "uid": bp_details.get("uid", ""),
                "title": bp_details.get("title", "") or s["title"],  # fallback to SRA title
                "description": bp_details.get("description", ""),
                "organism": bp_details.get("organism", ""),
                "data_type": bp_details.get("data_type", ""),
                "submission_date": bp_details.get("submission_date", ""),
                "last_update": bp_details.get("last_update", ""),
                "center_name": bp_details.get("center_name", ""),
                "linked_pubmed_ids": bp_details.get("linked_pubmed_ids", []),
                "url": bp_details.get("ncbi", {}).get("bioproject_url", ""),
            },

            # A “handle” back to how we found it
            "representative_sra": {
                "uid": uid,
                "title": s["title"],
                "url": f"https://www.ncbi.nlm.nih.gov/sra/?term={uid}",
            },

            # Aggregated quick facts (useful for cards + faceting)
            "run_summary": runinfo,

            # Provenance for append-only tracking
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

    # records.json: full append-only record list (what your UI can render)
    with open(f"{DB_DIR}/records.json", "w", encoding="utf-8") as f:
        json.dump(all_records, f, ensure_ascii=False, indent=2)

    # bioprojects.json: a convenient keyed map for fast lookup in the UI
    bp_map: Dict[str, Any] = {}
    for r in all_records:
        bp = (r.get("bioproject") or {}).get("accession") or (r.get("bioproject", {}).get("accession") if isinstance(r.get("bioproject"), dict) else None)
        if not bp and isinstance(r.get("bioproject"), dict):
            bp = r["bioproject"].get("accession")
        if not bp and isinstance(r.get("bioproject"), str):
            bp = r.get("bioproject")
        # In the new schema bp is nested dict; handle robustly:
        if isinstance(r.get("bioproject"), dict):
            bp = r["bioproject"].get("accession") or bp
            if bp:
                bp_map[bp] = r["bioproject"]

    with open(f"{DB_DIR}/bioprojects.json", "w", encoding="utf-8") as f:
        json.dump(bp_map, f, ensure_ascii=False, indent=2)

    # index.json: small metadata file your UI can fetch first
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

    # Load caches once; write back at end (fast + fewer writes)
    bp_cache = _load_bioproject_cache()
    uid_cache = _load_bioproject_uid_cache()

    try:
        if args.cmd == "backfill-year":
            year = args.year
            # NOTE: naive 365-day loop (keeps your original behavior).
            # If you care about leap years, swap to a proper year-day iterator.
            for day in (dt.date(year, 1, 1) + dt.timedelta(n) for n in range(365)):
                ds = day.isoformat()
                added, ids, projs = run_day(
                    ds, args.query, args.max_per_day,
                    seen_ids, seen_projects,
                    bp_cache, uid_cache
                )
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
                added, ids, projs = run_day(
                    ds, args.query, args.max_per_day,
                    seen_ids, seen_projects,
                    bp_cache, uid_cache
                )
                if added:
                    append_jsonl(f"{DATA_DIR}/catalog_{end.year}.jsonl", added)
                append_lines(SEEN_IDS, ids)
                append_lines(SEEN_PROJECTS, projs)
                seen_ids |= ids
                seen_projects |= projs
                print(ds, len(added))

        rebuild_exports()

    finally:
        # Persist caches even if a day fails mid-run (helps long backfills)
        _save_bioproject_cache(bp_cache)
        _save_bioproject_uid_cache(uid_cache)

if __name__ == "__main__":
    main()
