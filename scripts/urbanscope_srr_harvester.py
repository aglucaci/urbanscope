#!/usr/bin/env python3
r"""
UrbanScope — Urban Microbiome SRA Harvester (SRR-flatten + Geo inference + Assay tagging + BioProject enrichment)

You asked for:
- One row per SRR ("all SRA files" mode) from SRA UIDs via RunInfo
- Geographic inference (country/city + lat/lon) from BioSample attributes (optional)
- Assay tagging (16S vs WGS vs RNA-seq vs Amplicon) from RunInfo + metadata text
- ALSO: fetch richer BioProject details (optional, cached) via E-utilities (no scraping)

Core design:
- Query NCBI SRA via E-utilities.
- For each NEW SRA UID, fetch RunInfo (CSV) and emit ONE RECORD PER SRR run.
- Optionally enrich each SRR:
    - BioSample XML attributes -> geo inference
    - BioProject esummary -> title/desc/center/dates/data_type + raw esummary items
- Append-only JSONL catalog per year + website-ready merged exports in docs/db/.

Outputs (SRR-first):
- data/srr_catalog_<YEAR>.jsonl
- docs/db/srr_records.json
- docs/db/srr_index.json
- docs/latest_srr.json
- (optional snapshots) docs/db/biosamples.json, docs/db/bioprojects.json
- debug files under data/debug/ and docs/debug/ when --debug

Author: Alexander G. Lucaci


# VIEW WEBSITE
python -m http.server 8000
http://localhost:8000/docs/

# DAILY CRAWl
python3 scripts/urbanscope_srr_harvester.py crawl --sort date --page-size 100 --stop-after-new-srr 200 --fetch-bioproject --fetch-biosample --debug
python3 scripts/urbanscope_srr_harvester.py crawl --sort date --page-size 100 --stop-after-new-srr 200 --fetch-bioproject --fetch-biosample


# BIG CRAWL 
python3 scripts\urbanscope_srr_harvester.py crawl --page-size 100 --max-total 50000 --fetch-biosample --fetch-bioproject --debug
python3 scripts\urbanscope_srr_harvester.py crawl --page-size 100 --max-total 50000 --fetch-biosample --fetch-bioproject



# DUMP 
python3 '.\scripts\urbanscope_srr_harvester.py' crawl --page-size 200 --max-total 100000 --fetch-biosample --fetch-bioproject
python3 scripts\urbanscope_srr_harvester.py crawl --page-size 100 --max-total 50000 --fetch-biosample --fetch-bioproject --debug
python3 scripts/urbanscope_srr_harvester.py daily --days 30 --max-per-day 500 --fetch-biosample --fetch-bioproject
python3 scripts/urbanscope_srr_harvester.py daily --recent-days 7 --days 7 --max-per-day 2000 --fetch-biosample --fetch-bioproject
python3 scripts/urbanscope_srr_harvester.py crawl \
  --sort date \
  --page-size 100 \
  --stop-after-new-srr 200 \
  --stop-after-empty-pages 10 \
  --fetch-bioproject --fetch-biosample --debug
python3 scripts/urbanscope_srr_harvester.py crawl --sort date --page-size 100 --stop-after-new-srr 200 --stop-after-empty-pages 10 --fetch-bioproject --fetch-biosample --debug



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
TOOL_NAME = os.getenv("NCBI_TOOL", "urbanscope-srr-harvester")
NCBI_EMAIL = os.getenv("NCBI_EMAIL", "")

DEFAULT_QUERY = (
    '((urban OR city OR cities OR subway OR transit OR "built environment" '
    'OR wastewater OR sewage OR stormwater OR "public transit" OR "surface swab") '
    'AND (microbiome OR metagenom* OR metatranscriptom* OR "shotgun metagenomic" '
    'OR "environmental swab" OR "RNA-seq" OR "amplicon"))'
)

DATA_DIR = "data"
DOCS_DIR = "docs"
DB_DIR = f"{DOCS_DIR}/db"
CACHE_DIR = f"{DATA_DIR}/cache"
DEBUG_DIR = f"{DATA_DIR}/debug"
DOCS_DEBUG_DIR = f"{DOCS_DIR}/debug"

SEEN_SRA_UIDS = f"{DATA_DIR}/seen_sra_uids.txt"   # prevents reprocessing a UID
SEEN_SRR_RUNS = f"{DATA_DIR}/seen_srr_runs.txt"   # prevents duplicate SRRs

# caches
BIOSAMPLE_CACHE = f"{CACHE_DIR}/biosample.json"            # BioSample accession -> parsed XML attrs
BIOPROJECT_CACHE = f"{CACHE_DIR}/bioproject.json"          # PRJ accession -> details
BIOPROJECT_UID_CACHE = f"{CACHE_DIR}/bioproject_uid.json"  # PRJ accession -> numeric UID

DOCS_LATEST_SRR = f"{DOCS_DIR}/latest_srr.json"
DOCS_LATEST_DEBUG = f"{DOCS_DEBUG_DIR}/latest_report.json"

BIOPROJECT_RE = re.compile(r"\bPRJ(?:NA|EB|DB)\d+\b", re.I)


# ------------------------
# Utilities
# ------------------------
def ensure_dirs():
    for p in [DATA_DIR, DOCS_DIR, DB_DIR, CACHE_DIR, DEBUG_DIR, DOCS_DEBUG_DIR]:
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
            with urllib.request.urlopen(req, timeout=60) as r:
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
    vals = [v for v in vals if v]
    if not vals:
        return
    with open(path, "a", encoding="utf-8") as f:
        for v in vals:
            f.write(v + "\n")

def read_json(path: str, default):
    if not os.path.exists(path):
        return default
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return default

def inc(d: Dict[str, int], k: str, n: int = 1):
    d[k] = d.get(k, 0) + n

def write_json(path: str, obj: Any):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)
    os.replace(tmp, path)

def append_jsonl(path: str, records: List[Dict[str, Any]]):
    if not records:
        return
    with open(path, "a", encoding="utf-8") as f:
        for r in records:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")

def append_jsonl_one(path: str, rec: Dict[str, Any]):
    with open(path, "a", encoding="utf-8") as f:
        f.write(json.dumps(rec, ensure_ascii=False) + "\n")

def iter_jsonl(path: str):
    if not os.path.exists(path):
        return
    with open(path, encoding="utf-8") as fh:
        for ln in fh:
            ln = ln.strip()
            if ln:
                yield json.loads(ln)

def _norm(s: str) -> str:
    return re.sub(r"\s+", " ", (s or "").strip().lower())


# ------------------------
# NCBI helpers
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

def esearch_day(db: str, term: str, day: str, retmax: int, datetype: str = "edat") -> Tuple[List[str], str]:
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
    ids = [x.text for x in root.findall(".//Id") if x.text]
    return ids, url

def esearch_history(db: str, term: str, retstart: int, retmax: int, sort: str = ""):
    params = _eutils_params({
        "db": db,
        "term": term,
        "retmode": "xml",
        "retstart": str(retstart),
        "retmax": str(retmax),
        "usehistory": "n",
    })
    if sort:
        params["sort"] = sort  # <-- THIS is the key line

    url = EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(params)
    root = parse_xml(http_get(url))

    ids = [x.text for x in root.findall(".//Id") if x.text]
    count_total = int((root.findtext(".//Count") or "0").strip() or "0")

    return ids, count_total, url

def esummary(db: str, ids: List[str]) -> Tuple[ET.Element, str]:
    if not ids:
        return ET.Element("EMPTY"), ""
    params = _eutils_params({"db": db, "id": ",".join(ids), "retmode": "xml"})
    url = EUTILS + "esummary.fcgi?" + urllib.parse.urlencode(params)
    return parse_xml(http_get(url)), url

def esearch_any(db: str, term: str, retmax: int = 10) -> Tuple[List[str], str]:
    params = _eutils_params({"db": db, "term": term, "retmode": "xml", "retmax": str(retmax)})
    url = EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(params)
    root = parse_xml(http_get(url))
    ids = [x.text for x in root.findall(".//Id") if x.text]
    return ids, url

def efetch_runinfo_text(uid: str) -> Tuple[str, str]:
    params = _eutils_params({"db": "sra", "id": uid, "rettype": "runinfo", "retmode": "text"})
    url = EUTILS + "efetch.fcgi?" + urllib.parse.urlencode(params)
    text = http_get(url).decode(errors="replace")
    return text, url

def esearch_recent(db: str, term: str, reldate_days: int, retmax: int, datetype: str = "edat") -> Tuple[List[str], str]:
    params = _eutils_params({
        "db": db,
        "term": term,
        "retmode": "xml",
        "reldate": str(reldate_days),
        "datetype": datetype,
        "retmax": str(retmax),
        "sort": "date",
    })
    url = EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(params)
    root = parse_xml(http_get(url))
    ids = [x.text for x in root.findall(".//Id") if x.text]
    return ids, url

def esummary_sra(uids: List[str]) -> Tuple[Dict[str, Dict[str, Any]], str]:
    """
    Pull Title + raw DocSum Items. Returns (uid -> dict, url_used).
    """
    if not uids:
        return {}, ""
    root, url = esummary("sra", uids)
    out: Dict[str, Dict[str, Any]] = {}
    for d in root.findall(".//DocSum"):
        uid = (d.findtext("Id") or "").strip()
        if not uid:
            continue
        items: Dict[str, Any] = {}
        for it in d.findall("Item"):
            name = it.attrib.get("Name", "")
            if not name:
                continue
            if list(it):
                sub = [x.text for x in it.findall(".//Item") if x.text]
                items[name] = sub if sub else (it.text or "").strip()
            else:
                items[name] = (it.text or "").strip()

        title = (items.get("Title") or "").strip()

        bioproject_guess = ""
        for v in items.values():
            if isinstance(v, str):
                m = BIOPROJECT_RE.search(v)
                if m:
                    bioproject_guess = m.group(0).upper()
                    break
            elif isinstance(v, list):
                for s in v:
                    m = BIOPROJECT_RE.search(s or "")
                    if m:
                        bioproject_guess = m.group(0).upper()
                        break

        out[uid] = {"uid": uid, "title": title, "bioproject_guess": bioproject_guess, "items": items}
    return out, url


# ------------------------
# RunInfo: flatten to SRR rows
# ------------------------
def parse_runinfo_rows(uid: str, max_rows: int = 200000) -> Tuple[List[Dict[str, str]], Dict[str, Any]]:
    text, url = efetch_runinfo_text(uid)
    rows = list(csv.DictReader(text.splitlines()))
    if max_rows and len(rows) > max_rows:
        rows = rows[:max_rows]
    cols = list(rows[0].keys()) if rows else []
    return rows, {"url": url, "columns": cols, "rows": len(rows)}


# ------------------------
# BioSample enrichment (geo + attributes)
# ------------------------
COUNTRY_HINTS = {
    "usa": "United States",
    "u.s.a": "United States",
    "united states": "United States",
    "uk": "United Kingdom",
    "u.k.": "United Kingdom",
    "england": "United Kingdom",
    "scotland": "United Kingdom",
    "uae": "United Arab Emirates",
}

def efetch_biosample_xml(accession_or_uid: str) -> Tuple[str, str]:
    params = _eutils_params({"db": "biosample", "id": accession_or_uid, "retmode": "xml"})
    url = EUTILS + "efetch.fcgi?" + urllib.parse.urlencode(params)
    xmltxt = http_get(url).decode(errors="replace")
    return xmltxt, url

def parse_biosample_attributes_from_xml(xmltxt: str) -> Dict[str, Any]:
    """
    Extract <Attribute attribute_name="...">value</Attribute> pairs.
    """
    out: Dict[str, Any] = {"attributes": {}}
    root = ET.fromstring(xmltxt.encode("utf-8", errors="ignore"))
    for attr in root.findall(".//Attribute"):
        key = (attr.attrib.get("attribute_name") or attr.attrib.get("harmonized_name") or "").strip()
        val = (attr.text or "").strip()
        if key and val:
            out["attributes"][key] = val
    title = root.findtext(".//Title") or ""
    organism = root.findtext(".//Organism/OrganismName") or root.findtext(".//OrganismName") or ""
    if title.strip():
        out["title"] = title.strip()
    if organism.strip():
        out["organism"] = organism.strip()
    return out

def get_biosample_details(biosample_accession: str, biosample_cache: Dict[str, Any]) -> Dict[str, Any]:
    biosample_accession = (biosample_accession or "").strip()
    if not biosample_accession:
        return {}
    if biosample_accession in biosample_cache:
        return biosample_cache[biosample_accession] or {}
    try:
        xmltxt, url = efetch_biosample_xml(biosample_accession)
        parsed = parse_biosample_attributes_from_xml(xmltxt)
        parsed["accession"] = biosample_accession
        parsed["efetch_url"] = url
        biosample_cache[biosample_accession] = parsed
        return parsed
    except Exception as e:
        biosample_cache[biosample_accession] = {"accession": biosample_accession, "error": str(e)}
        return biosample_cache[biosample_accession]

def infer_geo(biosample_details: Dict[str, Any], fallbacks: List[str]) -> Dict[str, Any]:
    attrs = (biosample_details or {}).get("attributes", {}) if isinstance(biosample_details, dict) else {}
    raw_geo = ""

    # common attribute keys
    for k in ["geo_loc_name", "geographic location", "geographic_location", "country", "location"]:
        if k in attrs and attrs[k]:
            raw_geo = str(attrs[k]).strip()
            break

    lat = lon = ""
    latlon = ""
    for k in ["lat_lon", "latitude and longitude", "latitude_longitude"]:
        if k in attrs and attrs[k]:
            latlon = str(attrs[k]).strip()
            break
    if latlon:
        m = re.search(r"(-?\d+(?:\.\d+)?)\s*[, ]\s*(-?\d+(?:\.\d+)?)", latlon)
        if m:
            lat, lon = m.group(1), m.group(2)

    country = city = region = ""
    if raw_geo:
        parts = [p.strip() for p in raw_geo.split(":")]
        if len(parts) >= 2:
            country = parts[0]
            rest = ":".join(parts[1:])
            bits = [b.strip() for b in rest.split(",") if b.strip()]
            if bits:
                city = bits[-1]
                if len(bits) >= 2:
                    region = bits[-2]
        else:
            bits = [b.strip() for b in raw_geo.split(",") if b.strip()]
            if bits:
                country = bits[0]
                if len(bits) >= 2:
                    city = bits[-1]
                if len(bits) >= 3:
                    region = bits[-2]

    c_norm = _norm(country)
    if c_norm in COUNTRY_HINTS:
        country = COUNTRY_HINTS[c_norm]
    elif country:
        country = " ".join(w.capitalize() for w in country.split())

    if not country and fallbacks:
        blob = _norm(" | ".join([x for x in fallbacks if x]))
        for k, v in COUNTRY_HINTS.items():
            if k in blob:
                country = v
                break

    return {
        "country": country,
        "region": region,
        "city": city,
        "lat": lat,
        "lon": lon,
        "raw": raw_geo,
        "biosample_accession": (biosample_details or {}).get("accession", ""),
    }


# ------------------------
# Assay tagging
# ------------------------
def classify_assay(run_row: Dict[str, str], sra_title: str, biosample_details: Dict[str, Any]) -> Dict[str, Any]:
    title = _norm(sra_title or "")
    attrs = (biosample_details or {}).get("attributes", {}) if isinstance(biosample_details, dict) else {}
    attr_blob = _norm(" ".join([f"{k}:{v}" for k, v in attrs.items()]))

    strat = _norm(run_row.get("LibraryStrategy", ""))
    src = _norm(run_row.get("LibrarySource", ""))
    sel = _norm(run_row.get("LibrarySelection", ""))

    blob = " | ".join([title, attr_blob, strat, src, sel])

    hits: List[str] = []
    tags: List[str] = []

    # Amplicon path
    if "amplicon" in strat or "amplicon" in blob:
        hits.append("amplicon")
        tags.append("amplicon")

        if "16s" in blob or ("rrna" in blob and "16s" in blob):
            hits.append("16s")
            tags.append("16S")
            return {"assay_class": "16S", "assay_tags": tags, "confidence": "high", "rationale": hits}
        if "its" in blob:
            hits.append("its")
            tags.append("ITS")
            return {"assay_class": "ITS", "assay_tags": tags, "confidence": "high", "rationale": hits}

        return {"assay_class": "Amplicon", "assay_tags": tags, "confidence": "high", "rationale": hits}

    # RNA-seq / metatranscriptome
    if strat in ("rna-seq", "transcriptome") or "rna-seq" in blob or "metatranscriptom" in blob:
        hits.append("rna-seq/metatranscriptome")
        tags.append("RNA")
        return {"assay_class": "RNA-seq", "assay_tags": tags, "confidence": "high", "rationale": hits}

    # WGS / shotgun / metagenomics
    if strat in ("wgs", "metagenomic") or "shotgun" in blob or "wgs" in blob or "metagenom" in blob:
        hits.append("wgs/shotgun/metagenomic")
        tags.append("shotgun")
        return {"assay_class": "WGS", "assay_tags": tags, "confidence": "high", "rationale": hits}

    # Weak targeted hints
    if "pcr" in sel or "rrna" in sel:
        hits.append("PCR/rRNA selection")
        tags.append("targeted")
        if "16s" in blob:
            hits.append("16s")
            tags.append("16S")
            return {"assay_class": "16S", "assay_tags": tags, "confidence": "medium", "rationale": hits}
        if "its" in blob:
            hits.append("its")
            tags.append("ITS")
            return {"assay_class": "ITS", "assay_tags": tags, "confidence": "medium", "rationale": hits}
        return {"assay_class": "Amplicon", "assay_tags": tags, "confidence": "medium", "rationale": hits}

    return {"assay_class": "Unknown", "assay_tags": [], "confidence": "low", "rationale": []}


# ------------------------
# BioProject enrichment (E-utilities, cached)
# ------------------------
def bioproject_accession_to_uid(accession: str, uid_cache: Dict[str, str]) -> Optional[str]:
    accession = (accession or "").strip().upper()
    if not accession:
        return None
    if accession in uid_cache:
        return uid_cache[accession] or None
    ids, _ = esearch_any("bioproject", f"{accession}[Accession]", retmax=5)
    uid = ids[0] if ids else None
    uid_cache[accession] = uid or ""
    return uid

def parse_bioproject_esummary(uid: str) -> Dict[str, Any]:
    """
    Pull rich BioProject fields via esummary + keep raw items.
    """
    root, _ = esummary("bioproject", [uid])
    docsum = root.find(".//DocSum")
    if docsum is None:
        return {"uid": uid}

    items: Dict[str, Any] = {}
    for it in docsum.findall("Item"):
        name = it.attrib.get("Name", "")
        if not name:
            continue
        if list(it):
            sub = [x.text for x in it.findall(".//Item") if x.text]
            items[name] = sub if sub else (it.text or "").strip()
        else:
            items[name] = (it.text or "").strip()

    acc = (items.get("Project_Acc") or items.get("Accession") or "").strip().upper()
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
            "bioproject_url": f"https://www.ncbi.nlm.nih.gov/bioproject/{uid}",
        },
        "esummary_items": items,
    }

def get_bioproject_details(accession: str,
                           bp_cache: Dict[str, Any],
                           uid_cache: Dict[str, str]) -> Dict[str, Any]:
    accession = (accession or "").strip().upper()
    if not accession:
        return {}
    if accession in bp_cache:
        return bp_cache[accession] or {}

    uid = bioproject_accession_to_uid(accession, uid_cache)
    if not uid:
        bp_cache[accession] = {"accession": accession, "uid": "", "error": "uid_not_found"}
        return bp_cache[accession]

    details = parse_bioproject_esummary(uid)
    details["accession"] = details.get("accession") or accession
    bp_cache[accession] = details
    return details


# ------------------------
# Debug paths
# ------------------------
def debug_paths(tag: str) -> Dict[str, str]:
    return {
        "report": f"{DEBUG_DIR}/report_{tag}.json",
        "decision": f"{DEBUG_DIR}/decision_log_{tag}.ndjson",
        "initial": f"{DEBUG_DIR}/initial_{tag}.json",
    }


# ------------------------
# Build SRR records for one SRA UID
# ------------------------
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

    out: List[Dict[str, Any]] = []
    for r in rows:
        srr = (r.get("Run") or "").strip()
        if not srr:
            continue

        # BioSample -> geo + extra text
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

        # BioProject enrichment (from RunInfo BioProject column; fallback to esummary guess)
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
            "runinfo_row": r,  # full per-run row = maximum metadata density
            "geo": geo,
            "assay": assay,
            "bioproject": bioproject_details,  # optional; denormalized convenience
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


# ------------------------
# Exports (SRR-first)
# ------------------------
def rebuild_srr_exports():
    years = sorted(
        int(f.split("_")[2].split(".")[0])
        for f in os.listdir(DATA_DIR)
        if f.startswith("srr_catalog_") and f.endswith(".jsonl")
    )

    all_rows: List[Dict[str, Any]] = []
    for y in years:
        all_rows.extend(list(iter_jsonl(f"{DATA_DIR}/srr_catalog_{y}.jsonl")))

    with open(f"{DB_DIR}/srr_records.json", "w", encoding="utf-8") as f:
        json.dump(all_rows, f, ensure_ascii=False, indent=2)

    with open(f"{DB_DIR}/srr_index.json", "w", encoding="utf-8") as f:
        json.dump(
            {"generated": utc_now(), "total_srr_records": len(all_rows), "years": years},
            f, ensure_ascii=False, indent=2
        )

    # cache snapshots into docs/db (optional but useful for UI joins)
    bp_cache = read_json(BIOPROJECT_CACHE, {})
    if bp_cache:
        with open(f"{DB_DIR}/bioprojects.json", "w", encoding="utf-8") as f:
            json.dump(bp_cache, f, ensure_ascii=False, indent=2)

    bs_cache = read_json(BIOSAMPLE_CACHE, {})
    if bs_cache:
        with open(f"{DB_DIR}/biosamples.json", "w", encoding="utf-8") as f:
            json.dump(bs_cache, f, ensure_ascii=False, indent=2)


# ------------------------
# Ingest a list of SRA UIDs -> SRR records
# ------------------------
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
    counters: Dict[str, int] = {}
    paths = debug_paths(tag)

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
            srr_rows, _runinfo_dbg = build_srr_records_for_sra_uid(
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
                append_jsonl_one(paths["decision"], {
                    "uid": uid,
                    "decision": "error",
                    "error": str(e),
                })

    if sra_mark_seen:
        append_lines(SEEN_SRA_UIDS, sra_mark_seen)
    if srr_mark_seen:
        append_lines(SEEN_SRR_RUNS, srr_mark_seen)

    report = {"tag": tag, "generated_utc": utc_now(), "counters": counters, "debug_files": paths if debug else {}}
    if debug:
        write_json(paths["report"], report)

    return added_srr, report


# ------------------------
# CLI
# ------------------------
def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    # daily
    d = sub.add_parser("daily")
    d.add_argument("--days", type=int, default=1)
    d.add_argument("--query", default=DEFAULT_QUERY)
    d.add_argument("--max-per-day", type=int, default=500)
    d.add_argument("--debug", action="store_true")
    d.add_argument("--fetch-biosample", action="store_true", help="Enrich SRRs with BioSample XML attributes for geo inference.")
    d.add_argument("--fetch-bioproject", action="store_true", help="Enrich SRRs with BioProject esummary details (cached).")
    d.add_argument("--runinfo-max-rows", type=int, default=200000, help="Cap runinfo rows per SRA UID (safety).")
    d.add_argument("--recent-days", type=int, default=7, help="Search SRA entries from the last N days (Entrez date).")
    # daily subparser (optional but useful)
    d.add_argument(
        "--sort",
        default="date",
        help="NCBI esearch sort order (use 'date' for newest-first)."
    )

    # backfill-year
    b = sub.add_parser("backfill-year")
    b.add_argument("--year", type=int, required=True)
    b.add_argument("--query", default=DEFAULT_QUERY)
    b.add_argument("--max-per-day", type=int, default=500)
    b.add_argument("--debug", action="store_true")
    b.add_argument("--fetch-biosample", action="store_true")
    b.add_argument("--fetch-bioproject", action="store_true")
    b.add_argument("--runinfo-max-rows", type=int, default=200000)

    # crawl
    c = sub.add_parser("crawl")
    c.add_argument("--query", default=DEFAULT_QUERY)
    c.add_argument("--page-size", type=int, default=500)
    c.add_argument("--max-total", type=int, default=0, help="0 = no cap (crawl all matching).")
    c.add_argument("--debug", action="store_true")
    c.add_argument("--fetch-biosample", action="store_true")
    c.add_argument("--fetch-bioproject", action="store_true")
    c.add_argument("--runinfo-max-rows", type=int, default=200000)
    c.add_argument("--stop-after-new-srr", type=int, default=0, help="Stop crawling after N new SRR records have been added (0 = disabled).")
    # crawl subparser
    c.add_argument(
        "--sort",
        default="date",
        help="NCBI esearch sort order (use 'date' for newest-first)."
    )

    args = ap.parse_args()
    ensure_dirs()

    seen_sra = load_set(SEEN_SRA_UIDS)
    seen_srr = load_set(SEEN_SRR_RUNS)

    biosample_cache = read_json(BIOSAMPLE_CACHE, {})
    bp_cache = read_json(BIOPROJECT_CACHE, {})
    bp_uid_cache = read_json(BIOPROJECT_UID_CACHE, {})

    latest_added: List[Dict[str, Any]] = []
    reports: List[Dict[str, Any]] = []

    try:
        if args.cmd == "daily":
            # One recent-window search (reliable) + dedupe handles "already seen"
            tag = f"recent_{args.recent_days}d"
            paths = debug_paths(tag)
        
            uids, esearch_url = esearch_recent(
                "sra",
                args.query,
                args.recent_days,
                args.max_per_day,
                datetype="edat",
            )
            summaries, esummary_url = esummary_sra(uids[:500])
        
            if args.debug:
                write_json(paths["initial"], {
                    "tag": tag,
                    "query": args.query,
                    "recent_days": args.recent_days,
                    "max_per_day": args.max_per_day,
                    "esearch_url": esearch_url,
                    "esummary_url": esummary_url,
                    "uids_count": len(uids),
                    "uids_sample": uids[:50],
                })
                append_jsonl_one(paths["decision"], {
                    "tag": tag,
                    "decision": "sra_esearch_recent",
                    "esearch_url": esearch_url,
                    "esummary_url": esummary_url,
                    "uids_count": len(uids),
                })
        
            if not uids:
                # Make it obvious what's happening
                print(tag, "NO UIDS returned.")
                print("esearch_url:", esearch_url)
                # Still write latest report so UI has something
                report = {
                    "tag": tag,
                    "generated_utc": utc_now(),
                    "counters": {"uids_input": 0, "uids_new": 0, "srr_emitted": 0},
                    "urls": {"esearch": esearch_url, "esummary": esummary_url},
                }
                reports.append(report)
            else:
                added_srr, report = ingest_uids_to_srr(
                    tag=tag,
                    uids=uids,
                    summaries=summaries,
                    biosample_cache=biosample_cache,
                    bp_cache=bp_cache,
                    bp_uid_cache=bp_uid_cache,
                    seen_sra=seen_sra,
                    seen_srr=seen_srr,
                    fetch_biosample=args.fetch_biosample,
                    fetch_bioproject=args.fetch_bioproject,
                    debug=args.debug,
                    runinfo_max_rows=args.runinfo_max_rows,
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
                print(tag,
                      f"uids={report['counters'].get('uids_input',0)}",
                      f"new_uids={report['counters'].get('uids_new',0)}",
                      f"srr_emitted={report['counters'].get('srr_emitted',0)}")
                if report.get("debug_files"):
                    print("debug initial:", report["debug_files"].get("initial", ""))

        elif args.cmd == "backfill-year":
            year = args.year
            start = dt.date(year, 1, 1)
            for n in range(366):
                day = start + dt.timedelta(n)
                if day.year != year:
                    break
                ds = day.isoformat()

                uids, esearch_url = esearch_day("sra", args.query, ds, args.max_per_day, datetype="edat")
                summaries, esummary_url = esummary_sra(uids[:500])

                if args.debug:
                    append_jsonl_one(debug_paths(ds)["decision"], {
                        "tag": ds,
                        "decision": "sra_esearch",
                        "esearch_url": esearch_url,
                        "esummary_url": esummary_url,
                        "uids_count": len(uids),
                    })

                added_srr, report = ingest_uids_to_srr(
                    tag=ds,
                    uids=uids,
                    summaries=summaries,
                    biosample_cache=biosample_cache,
                    bp_cache=bp_cache,
                    bp_uid_cache=bp_uid_cache,
                    seen_sra=seen_sra,
                    seen_srr=seen_srr,
                    fetch_biosample=args.fetch_biosample,
                    fetch_bioproject=args.fetch_bioproject,
                    debug=args.debug,
                    runinfo_max_rows=args.runinfo_max_rows,
                )

                if added_srr:
                    append_jsonl(f"{DATA_DIR}/srr_catalog_{year}.jsonl", added_srr)

                reports.append(report)
                print(ds,
                      f"uids={report['counters'].get('uids_input',0)}",
                      f"new_uids={report['counters'].get('uids_new',0)}",
                      f"srr_emitted={report['counters'].get('srr_emitted',0)}")

        else:  # crawl
            retstart = 0
            page = 0
            total_seen = 0
            max_total = int(args.max_total or 0)
            new_srr_total = 0
            
            while True:
                ids, count_total, esearch_url = esearch_history(
                        "sra",
                        args.query,
                        retstart=retstart,
                        retmax=args.page_size,
                        sort=args.sort,   # <-- NEW
                    )
                if not ids:
                    break

                summaries, esummary_url = esummary_sra(ids[:500])

                tag = f"crawl_{page:06d}"
                if args.debug:
                    write_json(debug_paths(tag)["initial"], {
                        "tag": tag,
                        "page": page,
                        "retstart": retstart,
                        "page_size": args.page_size,
                        "count_total": count_total,
                        "esearch_url": esearch_url,
                        "esummary_url": esummary_url,
                        "ids_count": len(ids),
                        "ids_sample": ids[:25],
                    })

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
                        debug=args.debug,
                        runinfo_max_rows=args.runinfo_max_rows,
                    )
                    
                    new_srr_count = report["counters"].get("srr_emitted", 0)
                    new_srr_total += new_srr_count
                    
                    if args.stop_after_new_srr and new_srr_total >= args.stop_after_new_srr:
                        print(f"Stopping crawl: reached {new_srr_total} new SRRs "
                              f"(limit={args.stop_after_new_srr})")
                        break

                    print(tag,
                        f"page_ids={len(ids)}",
                        f"count_total={count_total}",
                        f"new_srr_page={new_srr_count}",
                        f"new_srr_total={new_srr_total}",
                        f"processed_uids={total_seen}",
                    )

                    this_year = dt.date.today().year
                
                    if added_srr:
                        append_jsonl(f"{DATA_DIR}/srr_catalog_{this_year}.jsonl", added_srr)
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

                    total_seen += len(ids)
                    print(tag,
                          f"page_ids={len(ids)}",
                          f"count_total={count_total}",
                          f"new_uids={report['counters'].get('uids_new',0)}",
                          f"srr_emitted={report['counters'].get('srr_emitted',0)}",
                          f"processed={total_seen}")
                    
                    print(tag,
                        f"sort={args.sort}",
                        f"page_ids={len(ids)}",
                        f"count_total={count_total}",
                    )

                    retstart += args.page_size
                    page += 1
                    if retstart >= count_total:
                        break
                    if max_total and total_seen >= max_total:
                        break

        # Latest + reports for UI
        write_json(DOCS_LATEST_SRR, {
            "generated_utc": utc_now(),
            "count": len(latest_added),
            "items": latest_added[:5000],
        })
        write_json(DOCS_LATEST_DEBUG, {
            "generated_utc": utc_now(),
            "reports": reports,
        })

        rebuild_srr_exports()

    finally:
        write_json(BIOSAMPLE_CACHE, biosample_cache)
        write_json(BIOPROJECT_CACHE, bp_cache)
        write_json(BIOPROJECT_UID_CACHE, bp_uid_cache)


if __name__ == "__main__":
    main()
