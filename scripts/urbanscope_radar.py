#!/usr/bin/env python3
"""
UrbanScope — Urban Metagenomics & Metatranscriptomics Dataset Radar (GitHub-only)

Adds:
- City / country extraction (from SRA SAMPLE attributes when present; heuristic fallback)
- Study-type classification: air / wastewater / surface (keyword + attribute-based)
- Archive browsing by year (docs/archive.json)
- CSV export for collaborators (data/catalog_YYYY.csv + docs/latest.csv)

Stores:
- data/seen_ids.txt                : dedupe (sra:* and pmid:*)
- data/catalog_<YEAR>.jsonl        : append-only yearly catalog
- data/catalog_<YEAR>.csv          : append-only yearly CSV for collaborators
- docs/latest.json                 : latest added records (small)
- docs/latest.csv                  : latest added records (CSV)
- docs/archive.json                : year index + raw links to yearly catalogs

Modes:
  backfill-year --year YYYY
  daily --days N
"""

from __future__ import annotations
import argparse
import csv
import datetime as dt
import json
import os
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from typing import Dict, List, Any, Set, Tuple, Optional

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

DEFAULT_QUERY = (
    '(urban OR city OR cities OR subway OR transit OR "built environment" OR wastewater OR sewage OR aerosol OR airborne) '
    'AND (metagenomics OR metatranscriptomics OR "shotgun metagenomic" OR "RNA-seq")'
)

DATA_DIR = "data"
DOCS_DIR = "docs"
SEEN_PATH = os.path.join(DATA_DIR, "seen_ids.txt")
LATEST_JSON_PATH = os.path.join(DOCS_DIR, "latest.json")
LATEST_CSV_PATH = os.path.join(DOCS_DIR, "latest.csv")
ARCHIVE_JSON_PATH = os.path.join(DOCS_DIR, "archive.json")


# -------------------------
# HTTP + NCBI helpers
# -------------------------
def http_get(url: str, timeout: int = 30) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "urbanscope/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def esearch(db: str, term: str, mindate: str, maxdate: str, retmax: int, sort: str) -> List[str]:
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
    params = {"dbfrom": "pubmed", "db": "sra", "id": ",".join(pmids), "retmode": "xml"}
    url = EUTILS + "elink.fcgi?" + urllib.parse.urlencode(params)
    root = ET.fromstring(http_get(url))
    out: Dict[str, List[str]] = {}
    for linkset in root.findall("LinkSet"):
        pmid = (linkset.findtext("IdList/Id") or "").strip()
        sra_ids = [x.text.strip() for x in linkset.findall(".//Link/Id") if x.text]
        if pmid and sra_ids:
            out[pmid] = sra_ids
    return out


def batched(lst: List[str], n: int) -> List[List[str]]:
    return [lst[i:i+n] for i in range(0, len(lst), n)]


def efetch_sra_xml(uids: List[str]) -> str:
    """
    Fetch richer SRA XML. This can include SAMPLE attributes like geo_loc_name.
    Keep batches small.
    """
    if not uids:
        return ""
    params = {"db": "sra", "id": ",".join(uids), "retmode": "xml"}
    url = EUTILS + "efetch.fcgi?" + urllib.parse.urlencode(params)
    return http_get(url).decode("utf-8", errors="replace")


# -------------------------
# Parsing + enrichment
# -------------------------
KEY_AIR = re.compile(r"\b(air|airborne|aerosol|aerosols|bioaerosol|pm2\.5|particulate)\b", re.I)
KEY_WW = re.compile(r"\b(wastewater|sewage|sewer|influent|effluent|sludge|wwtp|treatment plant)\b", re.I)
KEY_SURF = re.compile(r"\b(surface|swab|fomite|touch|handrail|rail|doorknob|bench|seat|turnstile)\b", re.I)

def classify_study_type(text_blob: str) -> str:
    """
    Priority: wastewater > air > surface (you can change priority).
    """
    t = text_blob or ""
    if KEY_WW.search(t):
        return "wastewater"
    if KEY_AIR.search(t):
        return "air"
    if KEY_SURF.search(t):
        return "surface"
    return "unknown"


def parse_geo(geo_loc_name: str) -> Tuple[str, str]:
    """
    Tries to infer (city, country) from common formats:
      - "Country: City"
      - "City, Country"
      - "City - Country"
    Returns ("", "") if unknown.
    """
    if not geo_loc_name:
        return "", ""

    s = geo_loc_name.strip()

    # Country: City (NCBI BioSample common)
    if ":" in s:
        left, right = s.split(":", 1)
        country = left.strip()
        city = right.strip()
        return city, country

    # City, Country
    if "," in s:
        parts = [p.strip() for p in s.split(",") if p.strip()]
        if len(parts) >= 2:
            city = parts[0]
            country = parts[-1]
            return city, country

    # City - Country
    if " - " in s:
        parts = [p.strip() for p in s.split(" - ") if p.strip()]
        if len(parts) >= 2:
            city = parts[0]
            country = parts[-1]
            return city, country

    # Last-resort: if it’s a single token, treat as unknown
    return "", ""


def extract_text(node: Optional[ET.Element]) -> str:
    return (node.text or "").strip() if node is not None else ""


def sra_uid_to_link(uid: str) -> str:
    return f"https://www.ncbi.nlm.nih.gov/sra/?term={uid}" if uid else ""


def parse_sra_records_from_efetch(xml_text: str) -> Dict[str, Dict[str, Any]]:
    """
    Parse a subset of fields from SRA efetch XML.
    Returns dict keyed by sra_uid (the UID we requested) with enriched fields.
    """
    out: Dict[str, Dict[str, Any]] = {}
    if not xml_text.strip():
        return out

    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return out

    # NCBI SRA efetch returns EXPERIMENT_PACKAGE_SET/EXPERIMENT_PACKAGE
    for pkg in root.findall(".//EXPERIMENT_PACKAGE"):
        # The efetch doc doesn't always expose the numeric UID as a simple field;
        # but our request uses numeric IDs. We can still capture useful text fields,
        # and later we attach it to the UID we requested by ordering (best-effort).
        # However, often there is <STUDY>/<IDENTIFIERS>/<PRIMARY_ID> like SRPxxxx,
        # and <SAMPLE>/<IDENTIFIERS>/<PRIMARY_ID> like SRSxxxx, and <RUN_SET>/<RUN>/<IDENTIFIERS>/<PRIMARY_ID> SRRxxxx.
        study_title = extract_text(pkg.find(".//STUDY/DESCRIPTOR/STUDY_TITLE"))
        study_abs = extract_text(pkg.find(".//STUDY/DESCRIPTOR/STUDY_ABSTRACT"))
        sample_title = extract_text(pkg.find(".//SAMPLE/TITLE"))
        library_strategy = extract_text(pkg.find(".//LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY"))
        library_source = extract_text(pkg.find(".//LIBRARY_DESCRIPTOR/LIBRARY_SOURCE"))
        platform = extract_text(pkg.find(".//PLATFORM/*/INSTRUMENT_MODEL"))

        # SAMPLE attributes (often contains geo_loc_name)
        attrs = {}
        for attr in pkg.findall(".//SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE"):
            tag = extract_text(attr.find("TAG")).lower()
            val = extract_text(attr.find("VALUE"))
            if tag and val:
                attrs[tag] = val

        geo_loc_name = attrs.get("geo_loc_name", "") or attrs.get("geographic location", "")
        city, country = parse_geo(geo_loc_name)

        # If city/country not found, try a couple other tags
        if not country:
            country = attrs.get("country", "") or attrs.get("country of origin", "")
        if not city:
            city = attrs.get("city", "") or attrs.get("collection_site", "")

        # Build a blob for classification
        blob = " ".join([
            study_title, study_abs, sample_title, library_strategy, library_source, platform,
            geo_loc_name, json.dumps(attrs, ensure_ascii=False)
        ])
        study_type = classify_study_type(blob)

        rec = {
            "study_title": study_title,
            "study_abstract": study_abs[:400] + ("…" if len(study_abs) > 400 else ""),
            "sample_title": sample_title,
            "library_strategy": library_strategy,
            "library_source": library_source,
            "platform_model": platform,
            "geo_loc_name": geo_loc_name,
            "city": city,
            "country": country,
            "study_type": study_type,
        }

        # We need to key this record — attempt to use an SRR primary id if present (not numeric UID),
        # but the caller ultimately maps by numeric UID. We'll just stash as "enrichment" and
        # attach per UID in best-effort order.
        # For stability, return a list later; here we just push into out under a synthetic key.
        out_key = f"pkg_{len(out)}"
        out[out_key] = rec

    return out


def enrich_uids(uids: List[str]) -> Dict[str, Dict[str, Any]]:
    """
    Enrich SRA numeric IDs with efetch data. Because efetch XML doesn't reliably echo numeric UID,
    we do best-effort attachment by batch order. This is not perfect, but it’s useful and improves over time.

    If you later want perfect mapping, we can switch to using SRP/SRS/SRR identifiers instead of numeric UID.
    """
    enriched: Dict[str, Dict[str, Any]] = {}
    if not uids:
        return enriched

    for batch in batched(uids, 20):
        xml_text = efetch_sra_xml(batch)
        parsed = parse_sra_records_from_efetch(xml_text)

        # Best-effort mapping: assign parsed packages in order to uids in this batch
        pkgs = list(parsed.values())
        for i, uid in enumerate(batch):
            if i < len(pkgs):
                enriched[uid] = pkgs[i]
            else:
                enriched[uid] = {}

        time.sleep(0.34)

    return enriched


# -------------------------
# Storage helpers
# -------------------------
def ensure_dirs() -> None:
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(DOCS_DIR, exist_ok=True)


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


def catalog_paths_for_year(year: int) -> Tuple[str, str]:
    return (
        os.path.join(DATA_DIR, f"catalog_{year}.jsonl"),
        os.path.join(DATA_DIR, f"catalog_{year}.csv"),
    )


def append_jsonl(path: str, records: List[Dict[str, Any]]) -> None:
    if not records:
        return
    with open(path, "a", encoding="utf-8") as f:
        for r in records:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")


CSV_FIELDS = [
    "ingested_utc",
    "day_utc",
    "query",
    "sra_uid",
    "sra_link",
    "pmid",
    "paper",
    "study_title",
    "sample_title",
    "library_strategy",
    "library_source",
    "platform_model",
    "geo_loc_name",
    "city",
    "country",
    "study_type",
    "title_fallback",
]


def append_csv(path: str, records: List[Dict[str, Any]]) -> None:
    if not records:
        return
    file_exists = os.path.exists(path)
    with open(path, "a", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        if not file_exists:
            w.writeheader()
        for r in records:
            w.writerow({k: r.get(k, "") for k in CSV_FIELDS})


def write_latest(outputs: List[Dict[str, Any]], query: str) -> None:
    payload = {
        "generated_at_utc": dt.datetime.utcnow().isoformat(timespec="seconds"),
        "query": query,
        "added_records": len(outputs),
        "records": outputs[-200:],  # keep small for Pages
    }
    with open(LATEST_JSON_PATH, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)

    # latest.csv for collaborators
    with open(LATEST_CSV_PATH, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        w.writeheader()
        for r in outputs[-200:]:
            w.writerow({k: r.get(k, "") for k in CSV_FIELDS})


def count_lines(path: str) -> int:
    if not os.path.exists(path):
        return 0
    with open(path, "r", encoding="utf-8") as f:
        return sum(1 for _ in f)


def build_archive_json(owner: str, repo: str) -> None:
    """
    Creates docs/archive.json listing yearly catalogs and raw links.
    (Pages cannot serve /data by default when Pages source is /docs, so we link to GitHub raw.)
    """
    years = []
    for fname in os.listdir(DATA_DIR):
        m = re.match(r"catalog_(\d{4})\.jsonl$", fname)
        if not m:
            continue
        y = int(m.group(1))
        jsonl_path = os.path.join(DATA_DIR, fname)
        csv_path = os.path.join(DATA_DIR, f"catalog_{y}.csv")
        years.append({
            "year": y,
            "jsonl_records": count_lines(jsonl_path),
            "csv_records": max(0, count_lines(csv_path) - 1) if os.path.exists(csv_path) else 0,
            "jsonl_raw": f"https://raw.githubusercontent.com/{owner}/{repo}/main/data/catalog_{y}.jsonl",
            "csv_raw": f"https://raw.githubusercontent.com/{owner}/{repo}/main/data/catalog_{y}.csv",
        })

    years.sort(key=lambda x: x["year"], reverse=True)

    payload = {
        "generated_at_utc": dt.datetime.utcnow().isoformat(timespec="seconds"),
        "years": years,
    }
    with open(ARCHIVE_JSON_PATH, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)


# -------------------------
# Core pipeline
# -------------------------
def run_day(query: str, day: dt.date, retmax: int, seen: Set[str]) -> Tuple[List[Dict[str, Any]], List[str]]:
    """
    Returns dataset-level records (one per new SRA UID) + new seen IDs.
    Includes PubMed->SRA linkage as optional fields on the dataset record.
    """
    day_str = day.strftime("%Y/%m/%d")

    # 1) SRA direct search
    try:
        sra_uids = esearch("sra", query, day_str, day_str, retmax=retmax, sort="Date")
    except Exception:
        sra_uids = []
    time.sleep(0.34)

    # 2) PubMed search + elink to SRA
    try:
        pmids = esearch("pubmed", query, day_str, day_str, retmax=retmax, sort="pub+date")
    except Exception:
        pmids = []
    time.sleep(0.34)

    pmids_new = [p for p in pmids if f"pmid:{p}" not in seen]
    pmid_to_sra = elink_pubmed_to_sra(pmids_new[:200]) if pmids_new else {}
    time.sleep(0.34)

    # Collect candidate SRA uids, and track linkage pmid -> uids
    linked_uids = []
    uid_to_pmid: Dict[str, str] = {}
    for pmid, uids in pmid_to_sra.items():
        for uid in uids:
            linked_uids.append(uid)
            # first pmid wins; good enough for initial version
            uid_to_pmid.setdefault(uid, pmid)

    all_candidate_uids = list(dict.fromkeys(sra_uids + linked_uids))  # preserve order, unique
    new_uids = [u for u in all_candidate_uids if f"sra:{u}" not in seen]

    # Enrich new UIDs
    enrich_map = enrich_uids(new_uids)

    outputs: List[Dict[str, Any]] = []
    newly_seen: List[str] = []

    for uid in new_uids:
        e = enrich_map.get(uid, {}) or {}
        pmid = uid_to_pmid.get(uid, "")

        # Fallback blob (if enrichment missing) for classification
        title_fallback = e.get("study_title") or e.get("sample_title") or ""
        if not title_fallback:
            title_fallback = ""  # keep explicit

        blob = " ".join([
            e.get("study_title", ""), e.get("sample_title", ""), e.get("geo_loc_name", ""),
            e.get("library_strategy", ""), e.get("library_source", ""),
            e.get("platform_model", ""), title_fallback
        ])
        study_type = e.get("study_type") or classify_study_type(blob)

        rec = {
            "ingested_utc": dt.datetime.utcnow().isoformat(timespec="seconds"),
            "day_utc": day.isoformat(),
            "query": query,
            "sra_uid": uid,
            "sra_link": sra_uid_to_link(uid),

            "pmid": pmid,
            "paper": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "",

            "study_title": e.get("study_title", ""),
            "sample_title": e.get("sample_title", ""),
            "library_strategy": e.get("library_strategy", ""),
            "library_source": e.get("library_source", ""),
            "platform_model": e.get("platform_model", ""),

            "geo_loc_name": e.get("geo_loc_name", ""),
            "city": e.get("city", ""),
            "country": e.get("country", ""),

            "study_type": study_type,
            "title_fallback": title_fallback,
        }

        outputs.append(rec)
        newly_seen.append(f"sra:{uid}")
        if pmid:
            newly_seen.append(f"pmid:{pmid}")

    # Dedupe ids within-run and against seen
    out_ids = []
    local = set()
    for x in newly_seen:
        if x not in seen and x not in local:
            out_ids.append(x)
            local.add(x)

    return outputs, out_ids


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_back = sub.add_parser("backfill-year", help="Backfill one year (run year-by-year in GitHub Actions)")
    p_back.add_argument("--year
