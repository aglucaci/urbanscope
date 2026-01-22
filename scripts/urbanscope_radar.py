#!/usr/bin/env python3
"""
UrbanScope — Urban Metagenomics & Metatranscriptomics Radar (Loose + Collapse + PubMed visibility)

What this version does (per day):
- "Loose" ingestion: accept ANY SRA datasets returned by your SRA search
  (no title requirement; no dependence on SRA esummary BioProject fields).
- Resolve BioProject robustly:
    1) elink(sra -> bioproject UID) + esummary(bioproject UID -> PRJ accession)
    2) fallback: efetch(runinfo) BioProject column (PRJ* accession)
- Collapse all SRA hits to BioProjects (group-by PRJ accession).
- For each BioProject group:
    - pick a representative SRA UID
    - store group sizes (how many SRA UIDs matched that day)
    - (optional) store a sample of those UIDs
- PubMed: show exactly what’s being used (the elink URL) + what was returned initially (first N PMIDs).
- Debug outputs so you can see:
    - initial SRA IDs
    - initial summaries
    - SRA->BioProject mapping coverage
    - which BioProjects were seen vs newly added
    - why some SRA UIDs cannot be mapped to a BioProject

Outputs:
- data/debug/report_<DAY>.json
- data/debug/decision_log_<DAY>.ndjson
- docs/debug/latest_report.json
- docs/latest.json (new BioProjects added in this run)
- docs/db/records.json + docs/db/index.json + docs/db/bioprojects.json

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
CACHE_DIR = f"{DATA_DIR}/cache"
DEBUG_DIR = f"{DATA_DIR}/debug"
DOCS_DEBUG_DIR = f"{DOCS_DIR}/debug"

ARCHIVE_DIR = f"{DOCS_DIR}/archive"
ARCHIVE_CSV_DIR = f"{ARCHIVE_DIR}/csv"

SEEN_IDS = f"{DATA_DIR}/seen_ids.txt"
SEEN_PROJECTS = f"{DATA_DIR}/seen_bioprojects.txt"

BIOPROJECT_CACHE = f"{CACHE_DIR}/bioproject.json"          # accession -> details dict
BIOPROJECT_UID_CACHE = f"{CACHE_DIR}/bioproject_uid.json"  # accession -> uid (numeric)

DOCS_LATEST = f"{DOCS_DIR}/latest.json"
DOCS_LATEST_DEBUG = f"{DOCS_DEBUG_DIR}/latest_report.json"

BIOPROJECT_RE = re.compile(r"\bPRJ(?:NA|EB|DB)\d+\b", re.I)

# ------------------------
# Utilities
# ------------------------
def ensure_dirs():
    for p in [DATA_DIR, DOCS_DIR, DB_DIR, CACHE_DIR, DEBUG_DIR, DOCS_DEBUG_DIR, ARCHIVE_DIR, ARCHIVE_CSV_DIR]:
        os.makedirs(p, exist_ok=True)

def utc_now():
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
    if not os.path.exists(path): return set()
    return set(x.strip() for x in open(path, encoding="utf-8") if x.strip())

def append_lines(path: str, vals: Iterable[str]):
    vals = list(vals)
    if not vals: return
    with open(path, "a", encoding="utf-8") as f:
        for v in vals:
            f.write(v + "\n")

def append_jsonl(path: str, records: List[Dict]):
    if not records: return
    with open(path, "a", encoding="utf-8") as f:
        for r in records:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")

def append_jsonl_one(path: str, rec: Dict):
    with open(path, "a", encoding="utf-8") as f:
        f.write(json.dumps(rec, ensure_ascii=False) + "\n")

def iter_jsonl(path: str):
    if not os.path.exists(path): return
    with open(path, encoding="utf-8") as fh:
        for ln in fh:
            if ln.strip():
                yield json.loads(ln)

def read_json(path: str, default):
    if not os.path.exists(path): return default
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return default

def write_json(path: str, obj: Any):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)
    os.replace(tmp, path)

def chunked(xs: List[str], n: int):
    for i in range(0, len(xs), n):
        yield xs[i:i+n]

def inc(d: Dict[str, int], k: str, n: int = 1):
    d[k] = d.get(k, 0) + n

# ------------------------
# NCBI helpers
# ------------------------
def _eutils_params(extra: Dict[str, str]) -> Dict[str, str]:
    p = dict(extra)
    if NCBI_API_KEY: p["api_key"] = NCBI_API_KEY
    if TOOL_NAME: p["tool"] = TOOL_NAME
    if NCBI_EMAIL: p["email"] = NCBI_EMAIL
    return p

def esearch(db: str, term: str, day: str, retmax: int, datetype: str = "edat") -> Tuple[List[str], str]:
    """
    Returns (ids, url_used)
    NOTE: for SRA, datetype=edat tends to be more consistent than pdat.
    """
    params = _eutils_params({
        "db": db, "term": term, "retmode": "xml",
        "mindate": day, "maxdate": day, "datetype": datetype,
        "retmax": str(retmax)
    })
    url = EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(params)
    root = parse_xml(http_get(url))
    ids = [x.text for x in root.findall(".//Id") if x.text]
    return ids, url

def esearch_any(db: str, term: str, retmax: int = 20) -> Tuple[List[str], str]:
    params = _eutils_params({"db": db, "term": term, "retmode": "xml", "retmax": str(retmax)})
    url = EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(params)
    root = parse_xml(http_get(url))
    ids = [x.text for x in root.findall(".//Id") if x.text]
    return ids, url

def esummary(db: str, ids: List[str]) -> Tuple[ET.Element, str]:
    if not ids:
        return ET.Element("EMPTY"), ""
    params = _eutils_params({"db": db, "id": ",".join(ids), "retmode": "xml"})
    url = EUTILS + "esummary.fcgi?" + urllib.parse.urlencode(params)
    return parse_xml(http_get(url)), url

def elink(dbfrom: str, db: str, ids: List[str], linkname: Optional[str] = None) -> Tuple[ET.Element, str]:
    if not ids:
        return ET.Element("EMPTY"), ""
    params: Dict[str, str] = {"dbfrom": dbfrom, "db": db, "id": ",".join(ids), "retmode": "xml"}
    if linkname:
        params["linkname"] = linkname
    params = _eutils_params(params)
    url = EUTILS + "elink.fcgi?" + urllib.parse.urlencode(params)
    return parse_xml(http_get(url)), url

def esummary_sra(uids: List[str]) -> Tuple[Dict[str, Dict], str]:
    """
    Best-effort SRA summary: title + (unreliable) PRJ accession.
    Returns (summaries, url_used).
    """
    if not uids:
        return {}, ""
    root, url = esummary("sra", uids)
    out: Dict[str, Dict] = {}
    for d in root.findall(".//DocSum"):
        uid = (d.findtext("Id") or "").strip()
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
    return out, url

def efetch_runinfo_text(uid: str) -> Tuple[str, str]:
    params = _eutils_params({"db": "sra", "id": uid, "rettype": "runinfo", "retmode": "text"})
    url = EUTILS + "efetch.fcgi?" + urllib.parse.urlencode(params)
    text = http_get(url).decode(errors="replace")
    return text, url

def efetch_runinfo(uid: str) -> Dict[str, Any]:
    text, _ = efetch_runinfo_text(uid)
    rows = list(csv.DictReader(text.splitlines()))

    biosamples = sorted({r.get("BioSample") for r in rows if r.get("BioSample")})
    platforms: Dict[str, int] = {}
    lib_strategies: Dict[str, int] = {}
    for r in rows:
        p = r.get("Platform")
        if p: platforms[p] = platforms.get(p, 0) + 1
        ls = r.get("LibraryStrategy")
        if ls: lib_strategies[ls] = lib_strategies.get(ls, 0) + 1

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
# Robust BioProject resolution (SRA -> BioProject)
# ------------------------
def sra_uids_to_bioproject_uids(sra_uids: List[str]) -> Tuple[Dict[str, str], str]:
    """
    Map SRA UID -> BioProject UID using elink. Returns ({sra_uid: bioproject_uid}, url_used_last).
    """
    out: Dict[str, str] = {}
    last_url = ""
    for chunk in chunked(sra_uids, 200):
        root, url = elink("sra", "bioproject", chunk)
        last_url = url
        for linkset in root.findall(".//LinkSet"):
            sra_id = (linkset.findtext("./IdList/Id") or "").strip()
            bp_id = (linkset.findtext(".//LinkSetDb/Link/Id") or "").strip()
            if sra_id and bp_id:
                out[sra_id] = bp_id
    return out, last_url

def bioproject_uid_to_accession(uid: str) -> str:
    uid = (uid or "").strip()
    if not uid:
        return ""
    root, _ = esummary("bioproject", [uid])
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

def bioproject_from_runinfo(sra_uid: str) -> Tuple[str, Dict[str, Any]]:
    """
    Fallback: parse PRJ accession from SRA runinfo CSV BioProject column.
    Returns (accession, debug_info).
    """
    text, url = efetch_runinfo_text(sra_uid)
    rows = list(csv.DictReader(text.splitlines()))
    for r in rows:
        m = BIOPROJECT_RE.search((r.get("BioProject") or ""))
        if m:
            return m.group(0).upper(), {"source": "runinfo", "url": url}
    return "", {"source": "runinfo_missing", "url": url}

# ------------------------
# BioProject enrichment (cached) + PubMed visibility
# ------------------------
def bioproject_accession_to_uid(accession: str, uid_cache: Dict[str, str]) -> Optional[str]:
    accession = accession.upper()
    if accession in uid_cache:
        return uid_cache[accession] or None
    ids, _ = esearch_any("bioproject", f"{accession}[Accession]", retmax=5)
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
    }

def fetch_bioproject_details(accession: str,
                             bp_cache: Dict[str, Any],
                             uid_cache: Dict[str, str],
                             pubmed_limit: int = 25) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    Returns (details, pubmed_debug)
    pubmed_debug includes:
      - pubmed_elink_url
      - pmids_initial (first N)
      - pmid_count
    """
    accession = accession.upper()
    if accession in bp_cache:
        # still return pubmed_debug if present; else empty
        cached = bp_cache[accession]
        return cached, cached.get("_pubmed_debug", {})

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

    if uid:
        root, _ = esummary("bioproject", [uid])
        docsum = root.find(".//DocSum")
        if docsum is not None:
            parsed = parse_bioproject_esummary_docsum(docsum)
            parsed["accession"] = parsed.get("accession") or accession
            details.update(parsed)

    # PubMed via elink (this is the “PubMed search” being performed)
    pubmed_debug = {"pubmed_elink_url": "", "pmids_initial": [], "pmid_count": 0}
    if uid:
        try:
            link_root, link_url = elink("bioproject", "pubmed", [uid])
            pmids = sorted({x.text for x in link_root.findall(".//LinkSetDb/Link/Id") if x.text})
            details["linked_pubmed_ids"] = pmids[:200]
            pubmed_debug = {
                "pubmed_elink_url": link_url,
                "pmids_initial": pmids[:pubmed_limit],
                "pmid_count": len(pmids),
            }
        except Exception:
            pubmed_debug = {"pubmed_elink_url": "", "pmids_initial": [], "pmid_count": 0}

    # Store pubmed debug in cache (keeps the “what search ran” visible)
    details["_pubmed_debug"] = pubmed_debug
    bp_cache[accession] = details
    return details, pubmed_debug

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
# Core logic (Loose + Collapse-to-BioProject + Debug)
# ------------------------
def debug_paths(day: str) -> Dict[str, str]:
    return {
        "report": f"{DEBUG_DIR}/report_{day}.json",
        "decision": f"{DEBUG_DIR}/decision_log_{day}.ndjson",
        "initial": f"{DEBUG_DIR}/initial_{day}.json",
    }

def run_day(day: str,
            query: str,
            retmax: int,
            seen_projects: Set[str],
            bp_cache: Dict[str, Any],
            uid_cache: Dict[str, str],
            debug: bool = False,
            pubmed_limit: int = 25,
            keep_sra_uid_sample: int = 50) -> Tuple[List[Dict], Set[str], Dict[str, Any]]:

    counters: Dict[str, int] = {}
    paths = debug_paths(day)

    # 1) SRA search (loose: we accept all returned IDs)
    sra_ids, esearch_url = esearch("sra", query, day, retmax, datetype="edat")
    inc(counters, "sra_ids", len(sra_ids))

    # 2) SRA summaries (best-effort only; but we save them as “initial returned”)
    summaries, esummary_url = esummary_sra(sra_ids[:500])  # cap summary load
    inc(counters, "summaries", len(summaries))

    # 3) Map SRA UID -> BioProject UID
    sra_to_bpuid, elink_url = sra_uids_to_bioproject_uids(list(summaries.keys()) or sra_ids[:500])
    inc(counters, "elink_mapped", len(sra_to_bpuid))

    # 4) Resolve BioProject UID -> PRJ accession (collapse key)
    bpuid_to_acc: Dict[str, str] = {}
    for bpuid in sorted(set(sra_to_bpuid.values())):
        acc = bioproject_uid_to_accession(bpuid)
        if acc:
            bpuid_to_acc[bpuid] = acc
    inc(counters, "bpuid_resolved_to_accession", len(bpuid_to_acc))

    # Save “what came back initially”
    if debug:
        write_json(paths["initial"], {
            "day": day,
            "query": query,
            "esearch_url": esearch_url,
            "esummary_url": esummary_url,
            "elink_sra_to_bioproject_url": elink_url,
            "sra_ids_initial_count": len(sra_ids),
            "sra_ids_initial_sample": sra_ids[:50],
            "summaries_initial_count": len(summaries),
            "summaries_initial_sample": list(summaries.values())[:10],
            "sra_to_bpuid_mapped_count": len(sra_to_bpuid),
            "sra_to_bpuid_sample": dict(list(sra_to_bpuid.items())[:10]),
            "bpuid_to_accession_count": len(bpuid_to_acc),
            "bpuid_to_accession_sample": dict(list(bpuid_to_acc.items())[:10]),
        })

    # 5) Collapse: group SRA UIDs into BioProjects
    bp_groups: Dict[str, List[str]] = {}   # accession -> [sra_uid,...]
    unmapped_sra: List[str] = []

    # iterate over the SRA UIDs we have summaries for; if summaries empty, fall back to raw IDs
    candidate_uids = list(summaries.keys()) if summaries else sra_ids
    for uid in candidate_uids:
        bp = ""
        method = ""
        # esummary best-effort
        if uid in summaries and summaries[uid].get("bioproject"):
            bp = summaries[uid]["bioproject"].upper()
            method = "esummary"
        # elink path
        if not bp:
            bpuid = sra_to_bpuid.get(uid, "")
            if bpuid:
                bp = bpuid_to_acc.get(bpuid, "").upper()
                method = "elink"
        # runinfo fallback
        runinfo_debug = {}
        if not bp:
            bp, runinfo_debug = bioproject_from_runinfo(uid)
            if bp:
                method = "runinfo"

        if not bp:
            unmapped_sra.append(uid)
            inc(counters, "drop_no_bioproject")
            if debug:
                append_jsonl_one(paths["decision"], {
                    "day": day, "uid": uid, "decision": "drop",
                    "reason": "no_bioproject", "bp_method": method or runinfo_debug.get("source", ""),
                    "runinfo_url": runinfo_debug.get("url","")
                })
            continue

        bp_groups.setdefault(bp, []).append(uid)
        inc(counters, f"bp_method_{method or 'unknown'}")

    inc(counters, "collapsed_bioprojects_total", len(bp_groups))
    inc(counters, "unmapped_sra_uids", len(unmapped_sra))

    # 6) Emit records for NEW BioProjects only (append-only), but still log SEEN ones.
    added_records: List[Dict] = []
    new_projects: Set[str] = set()

    for bp, uids in sorted(bp_groups.items(), key=lambda kv: len(kv[1]), reverse=True):
        if bp in seen_projects:
            inc(counters, "seen_bioproject_groups")
            if debug:
                append_jsonl_one(paths["decision"], {
                    "day": day, "decision": "seen",
                    "bp": bp, "group_size": len(uids),
                    "representative_uid": uids[0],
                })
            continue

        # NEW bioproject: enrich + add record
        rep_uid = uids[0]
        runinfo = efetch_runinfo(rep_uid)
        bp_details, pubmed_dbg = fetch_bioproject_details(bp, bp_cache, uid_cache, pubmed_limit=pubmed_limit)

        rec = {
            "bioproject": {
                "accession": bp,
                "uid": bp_details.get("uid", ""),
                "title": bp_details.get("title", "") or (summaries.get(rep_uid, {}).get("title","") if summaries else ""),
                "description": bp_details.get("description", ""),
                "organism": bp_details.get("organism", ""),
                "data_type": bp_details.get("data_type", ""),
                "submission_date": bp_details.get("submission_date", ""),
                "last_update": bp_details.get("last_update", ""),
                "center_name": bp_details.get("center_name", ""),
                "linked_pubmed_ids": bp_details.get("linked_pubmed_ids", []),
                "url": bp_details.get("ncbi", {}).get("bioproject_url", ""),
                "pubmed_debug": pubmed_dbg,   # <-- "what pubmed search is" + "what returned initially"
            },
            "collapsed": {
                "group_size_sra_uids": len(uids),
                "sra_uid_sample": uids[:keep_sra_uid_sample],
            },
            "representative_sra": {
                "uid": rep_uid,
                "title": summaries.get(rep_uid, {}).get("title", "") if summaries else "",
                "url": f"https://www.ncbi.nlm.nih.gov/sra/?term={rep_uid}",
            },
            "run_summary": runinfo,
            "provenance": {"source": "sra", "day": day, "ingested_utc": utc_now(), "query": query},
        }

        added_records.append(rec)
        new_projects.add(bp)
        inc(counters, "added_bioprojects")

        if debug:
            append_jsonl_one(paths["decision"], {
                "day": day, "decision": "add",
                "bp": bp, "group_size": len(uids),
                "representative_uid": rep_uid,
                "pubmed_elink_url": pubmed_dbg.get("pubmed_elink_url",""),
                "pmid_count": pubmed_dbg.get("pmid_count", 0),
                "pmids_initial": pubmed_dbg.get("pmids_initial", []),
            })

    report = {
        "day": day,
        "generated_utc": utc_now(),
        "query": query,
        "retmax": retmax,
        "urls": {
            "sra_esearch": esearch_url,
            "sra_esummary": esummary_url,
            "sra_to_bioproject_elink": elink_url,
        },
        "counters": counters,
        "top_bioproject_groups_sample": [
            {"bioproject": bp, "group_size": len(uids), "representative_uid": uids[0]}
            for bp, uids in list(sorted(bp_groups.items(), key=lambda kv: len(kv[1]), reverse=True))[:20]
        ],
        "unmapped_sra_uid_sample": unmapped_sra[:25],
        "debug_files": paths if debug else {},
    }

    if debug:
        write_json(paths["report"], report)

    return added_records, new_projects, report

# ------------------------
# Exports + latest files
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
    b.add_argument("--debug", action="store_true")
    b.add_argument("--pubmed-limit", type=int, default=25)

    d = sub.add_parser("daily")
    d.add_argument("--days", type=int, default=1)
    d.add_argument("--query", default=DEFAULT_QUERY)
    d.add_argument("--max-per-day", type=int, default=200)
    d.add_argument("--debug", action="store_true")
    d.add_argument("--pubmed-limit", type=int, default=25)

    args = ap.parse_args()
    ensure_dirs()

    seen_projects = load_set(SEEN_PROJECTS)
    seen_ids = load_set(SEEN_IDS)

    bp_cache = read_json(BIOPROJECT_CACHE, {})
    uid_cache = read_json(BIOPROJECT_UID_CACHE, {})

    latest_added: List[Dict[str, Any]] = []
    reports: List[Dict[str, Any]] = []

    try:
        if args.cmd == "backfill-year":
            year = args.year
            for day in (dt.date(year, 1, 1) + dt.timedelta(n) for n in range(365)):
                ds = day.isoformat()
                added, new_projects, report = run_day(
                    ds, args.query, args.max_per_day,
                    seen_projects, bp_cache, uid_cache,
                    debug=args.debug, pubmed_limit=args.pubmed_limit
                )

                if added:
                    append_jsonl(f"{DATA_DIR}/catalog_{year}.jsonl", added)

                # Update seen
                if new_projects:
                    append_lines(SEEN_PROJECTS, sorted(new_projects))
                    seen_projects |= new_projects

                # We still mark that we processed “day”
                reports.append(report)
                print(ds, f"added_bioprojects={report['counters'].get('added_bioprojects',0)}",
                      f"sra_ids={report['counters'].get('sra_ids',0)}",
                      f"collapsed_total={report['counters'].get('collapsed_bioprojects_total',0)}")

        else:
            end = dt.date.today()
            for i in range(args.days):
                ds = (end - dt.timedelta(i)).isoformat()
                added, new_projects, report = run_day(
                    ds, args.query, args.max_per_day,
                    seen_projects, bp_cache, uid_cache,
                    debug=args.debug, pubmed_limit=args.pubmed_limit
                )

                if added:
                    append_jsonl(f"{DATA_DIR}/catalog_{end.year}.jsonl", added)
                    for r in added:
                        bp = r.get("bioproject", {})
                        latest_added.append({
                            "day": ds,
                            "accession": bp.get("accession",""),
                            "title": bp.get("title",""),
                            "organism": bp.get("organism",""),
                            "data_type": bp.get("data_type",""),
                            "center_name": bp.get("center_name",""),
                            "url": bp.get("url",""),
                            "pubmed_elink_url": (bp.get("pubmed_debug") or {}).get("pubmed_elink_url",""),
                            "pmid_count": (bp.get("pubmed_debug") or {}).get("pmid_count", 0),
                            "pmids_initial": (bp.get("pubmed_debug") or {}).get("pmids_initial", []),
                        })

                if new_projects:
                    append_lines(SEEN_PROJECTS, sorted(new_projects))
                    seen_projects |= new_projects

                reports.append(report)
                print(ds, f"added_bioprojects={report['counters'].get('added_bioprojects',0)}",
                      f"sra_ids={report['counters'].get('sra_ids',0)}",
                      f"collapsed_total={report['counters'].get('collapsed_bioprojects_total',0)}")

        # Always write latest files so UI/debug has something to read
        write_json(DOCS_LATEST, {
            "generated_utc": utc_now(),
            "count": len(latest_added),
            "items": latest_added[:500],
        })
        write_json(DOCS_LATEST_DEBUG, {
            "generated_utc": utc_now(),
            "reports": reports,
        })

        rebuild_exports()

    finally:
        write_json(BIOPROJECT_CACHE, bp_cache)
        write_json(BIOPROJECT_UID_CACHE, uid_cache)

if __name__ == "__main__":
    main()
