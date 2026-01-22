#!/usr/bin/env python3
"""
UrbanScope — Urban Metagenomics & Metatranscriptomics Dataset Radar (GitHub-only)

You said:
- "gather all of this information and store it, I will parse it later in the html"
- "only show unique bioprojects"
- "include the title"
- city/country strict filter was too strict → make optional

This script does exactly that.

Core behavior
-------------
- Enforces **one record per BioProject** (PRJNA/PRJEB/PRJDB) across the entire database.
  - Records WITHOUT a BioProject accession are skipped.
- Always includes a **title** (from SRA where available; PubMed title also stored when linked).
- Stores rich metadata for later HTML parsing:
  - SRA: selected summary fields + a compact 'sra_items' dict (safe subset)
  - PubMed: title, journal, year, authors (first few), doi (if present)
  - Derived: study_type, city, country (best-effort), links
- Backfill year-by-year + daily incremental
- Writes GitHub Pages artifacts:
  - docs/latest.json (latest added bioproject records)
  - docs/db/index.json, docs/db/year/<YEAR>.json
  - docs/db/records.ndjson + docs/db/records.csv
  - docs/archive/* + CSV exports

Usage
-----
  # Backfill one year
  python scripts/urbanscope_radar.py backfill-year --year 2018 --max-per-day 200

  # Daily incremental
  python scripts/urbanscope_radar.py daily --days 1 --max-per-day 200

  # Optional strict location requirement (keep only records with BOTH city and country)
  python scripts/urbanscope_radar.py daily --days 1 --require-location

Notes
-----
- 100% stdlib
- NCBI E-utilities can disconnect; robust retry/backoff included.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import http.client
import json
import os
import random
import re
import socket
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple


EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

DEFAULT_QUERY = (
    '(urban OR city OR cities OR metropolis OR megacity OR subway OR transit OR "public transit" '
    'OR "built environment" OR "indoor microbiome" OR wastewater OR sewage OR stormwater) '
    'AND (metagenomics OR metatranscriptomics OR "shotgun metagenomic" OR "RNA-seq" OR "total RNA")'
)

# -------------------------
# Paths / folders
# -------------------------
DATA_DIR = "data"
DOCS_DIR = "docs"

ARCHIVE_DIR = os.path.join(DOCS_DIR, "archive")
ARCHIVE_CSV_DIR = os.path.join(ARCHIVE_DIR, "csv")

DB_DIR = os.path.join(DOCS_DIR, "db")
DB_YEAR_DIR = os.path.join(DB_DIR, "year")
DB_INDEX_JSON = os.path.join(DB_DIR, "index.json")
DB_ALL_NDJSON = os.path.join(DB_DIR, "records.ndjson")
DB_ALL_CSV = os.path.join(DB_DIR, "records.csv")

SEEN_IDS_PATH = os.path.join(DATA_DIR, "seen_ids.txt")                 # pmid:* and sra:*
SEEN_BIOPROJECTS_PATH = os.path.join(DATA_DIR, "seen_bioprojects.txt") # PRJ* only (uniqueness)

LATEST_JSON_PATH = os.path.join(DOCS_DIR, "latest.json")
DOCS_INDEX_HTML = os.path.join(DOCS_DIR, "index.html")

ARCHIVE_INDEX_JSON = os.path.join(ARCHIVE_DIR, "index.json")
ARCHIVE_INDEX_HTML = os.path.join(ARCHIVE_DIR, "index.html")
LATEST_CSV_PATH = os.path.join(ARCHIVE_CSV_DIR, "latest_added.csv")

# -------------------------
# Heuristics: study type, city/country
# -------------------------
STUDY_TYPE_RULES = [
    ("wastewater", [r"\bwastewater\b", r"\bsewage\b", r"\bstormwater\b", r"\binfluent\b", r"\beffluent\b", r"\bWWTP\b"]),
    ("air", [r"\bair\b", r"\bairborne\b", r"\baerosol\b", r"\baerosols\b", r"\bHVAC\b", r"\bventilation\b", r"\bdust\b"]),
    ("surface", [r"\bsurface\b", r"\bfomite\b", r"\btouch\b", r"\bhandrail\b", r"\bswab\b", r"\bdoorknob\b", r"\bbench\b"]),
]

COUNTRIES = sorted(
    {
        "United States", "USA", "United Kingdom", "UK", "England", "Scotland", "Wales", "Ireland",
        "Canada", "Mexico", "Brazil", "Argentina", "Chile", "Colombia", "Peru",
        "France", "Germany", "Italy", "Spain", "Portugal", "Netherlands", "Belgium", "Switzerland",
        "Austria", "Sweden", "Norway", "Denmark", "Finland", "Poland", "Czech Republic",
        "Greece", "Turkey", "Russia", "Ukraine",
        "China", "Hong Kong", "Taiwan", "Japan", "South Korea", "Korea", "Singapore",
        "India", "Pakistan", "Bangladesh", "Sri Lanka", "Nepal",
        "Thailand", "Vietnam", "Malaysia", "Indonesia", "Philippines",
        "Australia", "New Zealand",
        "South Africa", "Nigeria", "Kenya", "Ethiopia", "Egypt", "Morocco", "Ghana",
        "Israel", "Saudi Arabia", "United Arab Emirates", "UAE", "Qatar", "Iran", "Iraq",
    },
    key=lambda s: (-len(s), s.lower()),
)

CITY_ALIASES = {
    "New York City": [r"\bNYC\b", r"\bNew York City\b", r"\bNew York\b"],
    "London": [r"\bLondon\b"],
    "Paris": [r"\bParis\b"],
    "Tokyo": [r"\bTokyo\b"],
    "Beijing": [r"\bBeijing\b"],
    "Shanghai": [r"\bShanghai\b"],
    "Hong Kong": [r"\bHong Kong\b"],
    "Singapore": [r"\bSingapore\b"],
    "San Francisco": [r"\bSan Francisco\b", r"\bSF\b"],
    "Los Angeles": [r"\bLos Angeles\b", r"\bLA\b"],
    "Chicago": [r"\bChicago\b"],
    "Boston": [r"\bBoston\b"],
    "Seattle": [r"\bSeattle\b"],
    "Miami": [r"\bMiami\b"],
    "Toronto": [r"\bToronto\b"],
    "Vancouver": [r"\bVancouver\b"],
    "Montreal": [r"\bMontreal\b"],
    "Mexico City": [r"\bMexico City\b"],
    "São Paulo": [r"\bSao Paulo\b", r"\bSão Paulo\b"],
    "Rio de Janeiro": [r"\bRio de Janeiro\b"],
    "Buenos Aires": [r"\bBuenos Aires\b"],
    "Sydney": [r"\bSydney\b"],
    "Melbourne": [r"\bMelbourne\b"],
    "Auckland": [r"\bAuckland\b"],
    "Cape Town": [r"\bCape Town\b"],
    "Johannesburg": [r"\bJohannesburg\b"],
    "Nairobi": [r"\bNairobi\b"],
}

BIOPROJECT_RE = re.compile(r"\bPRJ(?:NA|EB|DB)\d+\b", re.IGNORECASE)


def utc_now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def classify_study_type(text: str) -> str:
    t = (text or "").lower()
    for label, patterns in STUDY_TYPE_RULES:
        for pat in patterns:
            if re.search(pat, t, flags=re.IGNORECASE):
                return label
    return "other"


def extract_city_country(text: str) -> Tuple[Optional[str], Optional[str]]:
    """Best-effort city/country extraction from title-like text (heuristic)."""
    if not text:
        return None, None

    country = None
    for c in COUNTRIES:
        if re.search(rf"\b{re.escape(c)}\b", text, flags=re.IGNORECASE):
            if c in {"USA", "United States"}:
                country = "United States"
            elif c in {"UK", "United Kingdom", "England", "Scotland", "Wales"}:
                country = "United Kingdom"
            elif c in {"UAE", "United Arab Emirates"}:
                country = "United Arab Emirates"
            elif c == "Korea":
                country = "South Korea"
            else:
                country = c
            break

    city = None
    for canonical, patterns in CITY_ALIASES.items():
        if any(re.search(pat, text, flags=re.IGNORECASE) for pat in patterns):
            city = canonical
            break

    # Heuristic: "City, Country"
    if not city:
        m = re.search(r"\b([A-Z][A-Za-z.\- ]{2,40}),\s*([A-Z][A-Za-z.\- ]{2,40})\b", text)
        if m:
            cand_city = m.group(1).strip()
            cand_country = m.group(2).strip()
            for c in COUNTRIES:
                if cand_country.lower() == c.lower() or re.search(rf"\b{re.escape(c)}\b", cand_country, flags=re.I):
                    if not country:
                        country = "United States" if c in {"USA", "United States"} else c
                    break
            if 2 <= len(cand_city) <= 40:
                city = cand_city

    return city, country


# -------------------------
# Robust HTTP for flaky upstream
# -------------------------
def http_get(url: str, timeout: int = 30, retries: int = 6, base_sleep: float = 0.8) -> bytes:
    """
    Robust GET with retries/backoff for NCBI E-utilities.
    Retries on RemoteDisconnected, timeouts, transient HTTP (429/5xx), URLError.
    """
    last_err: Optional[BaseException] = None

    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url,
                headers={
                    "User-Agent": "urbanscope/1.0 (GitHub Actions)",
                    "Accept": "application/xml,text/xml,application/json;q=0.9,*/*;q=0.8",
                },
            )
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = resp.read()
                if not data:
                    raise urllib.error.URLError("Empty response body")
                return data

        except (http.client.RemoteDisconnected, TimeoutError, socket.timeout, ConnectionResetError) as e:
            last_err = e
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504):
                last_err = e
            else:
                raise
        except urllib.error.URLError as e:
            last_err = e

        sleep_s = base_sleep * (2 ** attempt) + random.uniform(0, 0.25)
        time.sleep(sleep_s)

    raise RuntimeError(f"HTTP GET failed after {retries} retries. Last error: {last_err}")


def parse_xml(xml_bytes: bytes, url: str) -> ET.Element:
    try:
        return ET.fromstring(xml_bytes)
    except ET.ParseError as e:
        raise RuntimeError(f"XML parse error for URL: {url}") from e


# -------------------------
# NCBI E-utilities wrappers
# -------------------------
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
    root = parse_xml(http_get(url), url)
    return [n.text for n in root.findall(".//IdList/Id") if n.text]


def elink_pubmed_to_sra(pmids: List[str]) -> Dict[str, List[str]]:
    if not pmids:
        return {}
    params = {"dbfrom": "pubmed", "db": "sra", "id": ",".join(pmids), "retmode": "xml"}
    url = EUTILS + "elink.fcgi?" + urllib.parse.urlencode(params)
    root = parse_xml(http_get(url), url)

    out: Dict[str, List[str]] = {}
    for linkset in root.findall("LinkSet"):
        pmid = (linkset.findtext("IdList/Id") or "").strip()
        sra_ids = [x.text.strip() for x in linkset.findall(".//Link/Id") if x.text]
        if pmid and sra_ids:
            out[pmid] = sra_ids
    return out


def efetch_pubmed(pmids: List[str]) -> Dict[str, Dict[str, Any]]:
    """
    Fetch minimal PubMed metadata for later parsing:
    title, journal, year/medline date, doi, authors (first few).

    Returns {pmid: {...}}. Safe, compact.
    """
    if not pmids:
        return {}
    params = {"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"}
    url = EUTILS + "efetch.fcgi?" + urllib.parse.urlencode(params)
    root = parse_xml(http_get(url), url)

    out: Dict[str, Dict[str, Any]] = {}
    for art in root.findall(".//PubmedArticle"):
        pmid = (art.findtext(".//PMID") or "").strip()
        if not pmid:
            continue

        title = (art.findtext(".//ArticleTitle") or "").strip()
        journal = (art.findtext(".//Journal/Title") or "").strip()
        year = (art.findtext(".//PubDate/Year") or "").strip()
        medline = (art.findtext(".//PubDate/MedlineDate") or "").strip()
        pubdate = year or medline

        # DOI (best effort)
        doi = ""
        for aid in art.findall(".//ArticleIdList/ArticleId"):
            if (aid.attrib.get("IdType") or "").lower() == "doi":
                doi = (aid.text or "").strip()
                break

        # Authors (first 8)
        authors = []
        for a in art.findall(".//AuthorList/Author")[:8]:
            last = (a.findtext("LastName") or "").strip()
            fore = (a.findtext("ForeName") or "").strip()
            if fore and last:
                authors.append(f"{fore} {last}")
            elif last:
                authors.append(last)

        out[pmid] = {
            "pmid": pmid,
            "title": title,
            "journal": journal,
            "pubdate": pubdate,
            "doi": doi,
            "authors": authors,
            "paper": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
        }

    return out


def esummary_sra(uids: List[str]) -> Dict[str, Dict[str, Any]]:
    """
    Pull SRA summary for a set of uids and return mapping uid -> summary dict.

    We store a **safe subset** of useful fields in 'sra_items' for later HTML parsing.
    NCBI esummary can include large XML fields; we intentionally *do not* store ExpXml.
    """
    if not uids:
        return {}

    params = {"db": "sra", "id": ",".join(uids), "retmode": "xml"}
    url = EUTILS + "esummary.fcgi?" + urllib.parse.urlencode(params)
    root = parse_xml(http_get(url), url)

    out: Dict[str, Dict[str, Any]] = {}

    # Selected fields to capture if present (compact + useful)
    keep_names = {
        "Title",
        "Platform",
        "Strategy",
        "Source",
        "Selection",
        "TaxId",
        "Organism",
        "CreateDate",
        "UpdateDate",
        "Bioproject",
        "BioProject",
        "BioProjectAccn",
        "BioProjectAccession",
        "Study",
        "RunCount",
        "Bases",
        "BioprojectId",
        "Project",
    }

    for docsum in root.findall(".//DocSum"):
        uid = (docsum.findtext("Id") or "").strip()
        if not uid:
            continue

        sra_items: Dict[str, str] = {}
        title = ""
        bioproject = ""

        for it in docsum.findall("Item"):
            name = (it.attrib.get("Name", "") or "").strip()
            text = (it.text or "").strip()

            if name == "Title":
                title = text

            if name in keep_names and text:
                sra_items[name] = text

            # Try multiple BioProject fields
            if name.lower() in {"bioproject", "bioprojectaccn", "bioprojectaccession", "bioprojectacc", "project"}:
                m = BIOPROJECT_RE.search(text)
                if m:
                    bioproject = m.group(0).upper()

        # fallback: mine PRJ* from title
        if not bioproject and title:
            m = BIOPROJECT_RE.search(title)
            if m:
                bioproject = m.group(0).upper()

        out[uid] = {
            "sra_uid": uid,
            "title": title,
            "bioproject": bioproject,
            "sra_link": f"https://www.ncbi.nlm.nih.gov/sra/?term={uid}",
            "sra_items": sra_items,  # compact, later-parsable metadata
        }

    return out


# -------------------------
# Storage helpers
# -------------------------
def ensure_dirs() -> None:
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(DOCS_DIR, exist_ok=True)
    os.makedirs(ARCHIVE_DIR, exist_ok=True)
    os.makedirs(ARCHIVE_CSV_DIR, exist_ok=True)
    os.makedirs(DB_DIR, exist_ok=True)
    os.makedirs(DB_YEAR_DIR, exist_ok=True)


def load_lineset(path: str) -> Set[str]:
    if not os.path.exists(path):
        return set()
    with open(path, "r", encoding="utf-8") as f:
        return set(line.strip() for line in f if line.strip())


def append_lines(path: str, values: Iterable[str]) -> None:
    vals = [v for v in values if v]
    if not vals:
        return
    with open(path, "a", encoding="utf-8") as f:
        for v in vals:
            f.write(v + "\n")


def append_jsonl(path: str, records: List[Dict[str, Any]]) -> None:
    if not records:
        return
    with open(path, "a", encoding="utf-8") as f:
        for r in records:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")


def catalog_path_for_year(year: int) -> str:
    return os.path.join(DATA_DIR, f"catalog_{year}.jsonl")


def iter_jsonl(path: str) -> Iterable[Dict[str, Any]]:
    if not os.path.exists(path):
        return
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                yield json.loads(line)
            except json.JSONDecodeError:
                continue


def list_catalog_years() -> List[int]:
    years = []
    if not os.path.exists(DATA_DIR):
        return years
    for fn in os.listdir(DATA_DIR):
        m = re.match(r"catalog_(\d{4})\.jsonl$", fn)
        if m:
            years.append(int(m.group(1)))
    return sorted(years)


def daterange(start: dt.date, end: dt.date) -> Iterable[dt.date]:
    cur = start
    while cur <= end:
        yield cur
        cur += dt.timedelta(days=1)


# -------------------------
# BioProject uniqueness helpers
# -------------------------
def bioproject_key(summary: Dict[str, Any]) -> str:
    """
    Unique key = BioProject accession only.
    If missing → record is skipped.
    """
    raw = (summary.get("bioproject", "") or "").strip().upper()
    if not raw:
        return ""
    m = BIOPROJECT_RE.search(raw)
    return m.group(0).upper() if m else ""


# -------------------------
# Core acquisition (returns unique BioProject records)
# -------------------------
def run_day(
    query: str,
    day: dt.date,
    retmax: int,
    seen_ids: Set[str],
    seen_bioprojects: Set[str],
    require_location: bool,
) -> Tuple[List[Dict[str, Any]], List[str], Set[str]]:
    """
    Returns:
      added_records (one per new BioProject),
      new_seen_ids (pmid:* and sra:*),
      new_bioprojects (PRJ*)
    """
    mindate = day.strftime("%Y/%m/%d")
    maxdate = day.strftime("%Y/%m/%d")

    added: List[Dict[str, Any]] = []
    new_seen: List[str] = []
    new_projects: Set[str] = set()

    def nap() -> None:
        time.sleep(0.6)

    # -----------------
    # A) SRA direct search
    # -----------------
    try:
        sra_uids = esearch("sra", query, mindate, maxdate, retmax=retmax, sort="Date")
    except Exception as e:
        print(f"[{day.isoformat()}] SRA esearch ERROR: {e}")
        sra_uids = []
    nap()

    sra_uids_new = [u for u in sra_uids if f"sra:{u}" not in seen_ids]
    if sra_uids_new:
        # Summarize and decide which BioProjects are new
        sra_map = esummary_sra(sra_uids_new[:400])
        for uid, s in sra_map.items():
            bp = bioproject_key(s)
            if not bp:
                continue
            if bp in seen_bioprojects or bp in new_projects:
                continue

            title = (s.get("title") or "").strip()
            if not title:
                continue  # you asked to include the title

            stype = classify_study_type(title)
            city, country = extract_city_country(title)

            if require_location and (not city or not country):
                continue

            new_projects.add(bp)
            added.append(
                {
                    "bioproject": bp,
                    "representative": {
                        "sra_uid": uid,
                        "sra_link": s.get("sra_link", ""),
                        "title": title,
                    },
                    "derived": {
                        "study_type": stype,
                        "city": city,
                        "country": country,
                    },
                    "provenance": {
                        "source": "sra_search",
                        "day_utc": day.isoformat(),
                        "query": query,
                        "ingested_utc": utc_now_iso(),
                    },
                    # compact, later-parsable metadata
                    "sra_items": s.get("sra_items", {}),
                    # placeholder for later linking
                    "pubmed": None,
                    "linked_sra_uids": [uid],  # for this bioproject (we may expand below in PubMed path)
                }
            )
            new_seen.append(f"sra:{uid}")

    nap()

    # -----------------
    # B) PubMed -> SRA (choose 1 representative per BioProject)
    # -----------------
    # If a BioProject was already added via SRA-search above, PubMed-linked info can still enrich it,
    # but we keep it simple: only add new BioProjects not yet seen or added today.
    try:
        pmids = esearch("pubmed", query, mindate, maxdate, retmax=retmax, sort="pub+date")
    except Exception as e:
        print(f"[{day.isoformat()}] PubMed esearch ERROR: {e}")
        pmids = []
    nap()

    pmids_new = [p for p in pmids if f"pmid:{p}" not in seen_ids]
    if pmids_new:
        try:
            pmid_to_sra = elink_pubmed_to_sra(pmids_new[:250])
        except Exception as e:
            print(f"[{day.isoformat()}] PubMed elink ERROR: {e}")
            pmid_to_sra = {}

        # Fetch PubMed metadata for later parsing
        # Keep it compact: title/journal/pubdate/doi/authors
        pubmed_meta = efetch_pubmed(list(pmid_to_sra.keys())[:250]) if pmid_to_sra else {}

        # Summarize all linked SRA uids, group by BioProject
        all_uids = sorted({uid for uids in pmid_to_sra.values() for uid in uids})
        sra_map = esummary_sra(all_uids[:600]) if all_uids else {}

        nap()

        # Build mapping: pmid -> {bioproject -> [uids]}
        pmid_bp_uids: Dict[str, Dict[str, List[str]]] = {}
        for pmid, uids in pmid_to_sra.items():
            bp_map: Dict[str, List[str]] = {}
            for uid in uids:
                s = sra_map.get(uid, {"bioproject": ""})
                bp = bioproject_key(s)
                if not bp:
                    continue
                bp_map.setdefault(bp, []).append(uid)
            if bp_map:
                pmid_bp_uids[pmid] = bp_map

        for pmid, bp_map in pmid_bp_uids.items():
            # Choose one bioproject from this paper that is new
            # (If paper links multiple BioProjects, we may add more than one,
            #  but still strictly unique per BioProject)
            for bp, uids in bp_map.items():
                if bp in seen_bioprojects or bp in new_projects:
                    continue

                # pick representative uid with a non-empty title if possible
                chosen_uid = None
                chosen_summary: Optional[Dict[str, Any]] = None
                for uid in uids:
                    s = sra_map.get(uid)
                    if not s:
                        continue
                    if (s.get("title") or "").strip():
                        chosen_uid = uid
                        chosen_summary = s
                        break
                if not chosen_uid:
                    continue
                assert chosen_summary is not None

                title = (chosen_summary.get("title") or "").strip()
                if not title:
                    continue  # you asked to include the title

                stype = classify_study_type(title)
                city, country = extract_city_country(title)

                if require_location and (not city or not country):
                    continue

                new_projects.add(bp)

                added.append(
                    {
                        "bioproject": bp,
                        "representative": {
                            "sra_uid": chosen_uid,
                            "sra_link": chosen_summary.get("sra_link", f"https://www.ncbi.nlm.nih.gov/sra/?term={chosen_uid}"),
                            "title": title,
                        },
                        "derived": {
                            "study_type": stype,
                            "city": city,
                            "country": country,
                        },
                        "provenance": {
                            "source": "pubmed_elink_sra",
                            "day_utc": day.isoformat(),
                            "query": query,
                            "ingested_utc": utc_now_iso(),
                        },
                        "sra_items": chosen_summary.get("sra_items", {}),
                        "pubmed": pubmed_meta.get(pmid, {"pmid": pmid, "paper": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"}),
                        "linked_sra_uids": uids[:200],  # keep bounded
                    }
                )

                new_seen.append(f"pmid:{pmid}")
                new_seen.append(f"sra:{chosen_uid}")

    # De-dupe IDs within-run and against seen_ids
    uniq_ids: List[str] = []
    local = set()
    for x in new_seen:
        if x not in seen_ids and x not in local:
            uniq_ids.append(x)
            local.add(x)

    return added, uniq_ids, new_projects


# -------------------------
# Public artifacts (DB exports, CSV, archive, dashboard stub)
# -------------------------
CSV_FIELDS = [
    "bioproject",
    "rep_sra_uid",
    "rep_title",
    "study_type",
    "city",
    "country",
    "pmid",
    "paper",
    "doi",
    "journal",
    "pubdate",
    "authors",
    "source",
    "day_utc",
    "ingested_utc",
    "rep_sra_link",
]

def write_latest_json(records: List[Dict[str, Any]], query: str, require_location: bool) -> None:
    payload = {
        "generated_at_utc": utc_now_iso(),
        "query": query,
        "require_location": bool(require_location),
        "added_bioprojects": len(records),
        "records": records[-250:],
    }
    with open(LATEST_JSON_PATH, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)


def record_to_csv_row(r: Dict[str, Any]) -> Dict[str, Any]:
    rep = r.get("representative") or {}
    der = r.get("derived") or {}
    prov = r.get("provenance") or {}
    pm = r.get("pubmed") or {}

    return {
        "bioproject": r.get("bioproject", ""),
        "rep_sra_uid": rep.get("sra_uid", ""),
        "rep_title": rep.get("title", ""),
        "study_type": der.get("study_type", ""),
        "city": der.get("city", ""),
        "country": der.get("country", ""),
        "pmid": pm.get("pmid", ""),
        "paper": pm.get("paper", ""),
        "doi": pm.get("doi", ""),
        "journal": pm.get("journal", ""),
        "pubdate": pm.get("pubdate", ""),
        "authors": " ; ".join(pm.get("authors", []) or []) if isinstance(pm.get("authors"), list) else (pm.get("authors") or ""),
        "source": prov.get("source", ""),
        "day_utc": prov.get("day_utc", ""),
        "ingested_utc": prov.get("ingested_utc", ""),
        "rep_sra_link": rep.get("sra_link", ""),
    }


def write_latest_csv(records: List[Dict[str, Any]]) -> None:
    with open(LATEST_CSV_PATH, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        w.writeheader()
        for r in records[-5000:]:
            w.writerow(record_to_csv_row(r))


def build_archive_index() -> Dict[str, Any]:
    years = []
    for y in list_catalog_years():
        path = catalog_path_for_year(y)
        count = 0
        with open(path, "r", encoding="utf-8") as f:
            for _ in f:
                count += 1
        years.append({"year": y, "records": count, "jsonl": f"../../data/catalog_{y}.jsonl"})
    payload = {"generated_at_utc": utc_now_iso(), "years": years}
    with open(ARCHIVE_INDEX_JSON, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)
    return payload


def write_archive_html(index_payload: Dict[str, Any]) -> None:
    years = index_payload.get("years", [])
    rows = []
    for y in years:
        rows.append(f"<tr><td><b>{y['year']}</b></td><td>{y['records']}</td></tr>")
    html = f"""<!doctype html>
<html><head>
<meta charset="utf-8"/><meta name="viewport" content="width=device-width,initial-scale=1"/>
<title>UrbanScope Archive</title>
<style>
body{{font-family:system-ui,Segoe UI,Roboto,Arial;max-width:980px;margin:40px auto;padding:0 16px;line-height:1.55}}
.muted{{color:#666}} table{{border-collapse:collapse;width:100%;margin-top:14px}}
th,td{{border-bottom:1px solid #eee;padding:10px;text-align:left}}
code{{background:#f2f2f2;padding:2px 6px;border-radius:6px}}
</style></head><body>
<h1>UrbanScope Archive</h1>
<p class="muted">Index updated: <code>{index_payload.get("generated_at_utc","")}</code></p>
<p><a href="../index.html">← Back</a> · <a href="../db/index.json">db/index.json</a></p>
<table><thead><tr><th>Year</th><th>BioProjects</th></tr></thead>
<tbody>
{''.join(rows) if rows else '<tr><td colspan="2" class="muted">No archive yet.</td></tr>'}
</tbody></table>
</body></html>
"""
    with open(ARCHIVE_INDEX_HTML, "w", encoding="utf-8") as f:
        f.write(html)


def write_database_exports() -> Dict[str, Any]:
    years = list_catalog_years()
    total = 0
    by_year: Dict[str, int] = {}
    by_type: Dict[str, int] = {}
    by_country: Dict[str, int] = {}

    # Full NDJSON
    with open(DB_ALL_NDJSON, "w", encoding="utf-8") as nd:
        for y in years:
            for r in iter_jsonl(catalog_path_for_year(y)):
                nd.write(json.dumps(r, ensure_ascii=False) + "\n")
                total += 1
                by_year[str(y)] = by_year.get(str(y), 0) + 1
                st = (r.get("derived") or {}).get("study_type") or "other"
                by_type[st] = by_type.get(st, 0) + 1
                c = (r.get("derived") or {}).get("country") or "Unknown"
                by_country[c] = by_country.get(c, 0) + 1

    # Full CSV
    with open(DB_ALL_CSV, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        w.writeheader()
        for y in years:
            for r in iter_jsonl(catalog_path_for_year(y)):
                w.writerow(record_to_csv_row(r))

    # Per-year JSON slices (fast UI)
    for y in years:
        out = os.path.join(DB_YEAR_DIR, f"{y}.json")
        items = list(iter_jsonl(catalog_path_for_year(y)))
        with open(out, "w", encoding="utf-8") as f:
            json.dump({"year": y, "records": items}, f, ensure_ascii=False)

    index = {
        "generated_at_utc": utc_now_iso(),
        "years": years,
        "total_records": total,
        "counts_by_year": by_year,
        "counts_by_type": dict(sorted(by_type.items(), key=lambda x: x[0])),
        "top_countries": sorted(by_country.items(), key=lambda x: x[1], reverse=True)[:25],
        "downloads": {"csv": "db/records.csv", "ndjson": "db/records.ndjson", "archive": "archive/index.html"},
    }
    with open(DB_INDEX_JSON, "w", encoding="utf-8") as f:
        json.dump(index, f, indent=2, ensure_ascii=False)

    return index


def write_dashboard_stub() -> None:
    """
    Minimal stub; you said you'll parse later in HTML.
    This just links to db + latest artifacts.
    """
    html = """<!doctype html>
<html><head>
<meta charset="utf-8"/><meta name="viewport" content="width=device-width,initial-scale=1"/>
<title>UrbanScope</title>
<style>
body{font-family:system-ui,Segoe UI,Roboto,Arial;max-width:980px;margin:40px auto;padding:0 16px;line-height:1.55}
.muted{color:#666}
a{text-decoration:none} a:hover{text-decoration:underline}
code{background:#f2f2f2;padding:2px 6px;border-radius:6px}
</style></head><body>
<h1>UrbanScope</h1>
<p class="muted">Urban metagenomics & metatranscriptomics BioProject database (updated daily).</p>
<ul>
  <li><a href="db/index.json">db/index.json</a></li>
  <li><a href="db/records.ndjson">db/records.ndjson</a> (full database)</li>
  <li><a href="db/records.csv">db/records.csv</a> (collaborator export)</li>
  <li><a href="latest.json">latest.json</a> (new BioProjects this run)</li>
  <li><a href="archive/index.html">archive/index.html</a></li>
  <li><a href="archive/csv/latest_added.csv">archive/csv/latest_added.csv</a></li>
</ul>
</body></html>
"""
    with open(DOCS_INDEX_HTML, "w", encoding="utf-8") as f:
        f.write(html)


# -------------------------
# Main
# -------------------------
def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_back = sub.add_parser("backfill-year", help="Backfill one year (GitHub Actions-friendly)")
    p_back.add_argument("--year", type=int, required=True)
    p_back.add_argument("--query", default=DEFAULT_QUERY)
    p_back.add_argument("--max-per-day", type=int, default=200)
    p_back.add_argument("--require-location", action="store_true",
                        help="Only keep records with BOTH city and country (default: off)")

    p_daily = sub.add_parser("daily", help="Daily incremental update (last N days)")
    p_daily.add_argument("--days", type=int, default=1)
    p_daily.add_argument("--query", default=DEFAULT_QUERY)
    p_daily.add_argument("--max-per-day", type=int, default=200)
    p_daily.add_argument("--require-location", action="store_true",
                         help="Only keep records with BOTH city and country (default: off)")

    args = ap.parse_args()
    ensure_dirs()

    seen_ids = load_lineset(SEEN_IDS_PATH)
    seen_bioprojects = load_lineset(SEEN_BIOPROJECTS_PATH)

    added_all: List[Dict[str, Any]] = []
    new_ids_all: List[str] = []
    new_bps_all: Set[str] = set()

    if args.cmd == "backfill-year":
        year = int(args.year)
        start = dt.date(year, 1, 1)
        end = dt.date(year, 12, 31)
        cat_path = catalog_path_for_year(year)

        for day in daterange(start, end):
            try:
                added, new_ids, new_bps = run_day(
                    args.query, day, args.max_per_day, seen_ids, seen_bioprojects, bool(args.require_location)
                )
            except Exception as e:
                print(f"[{day.isoformat()}] ERROR: {e}")
                continue

            if added:
                append_jsonl(cat_path, added)
                added_all.extend(added)

            if new_ids:
                append_lines(SEEN_IDS_PATH, new_ids)
                for x in new_ids:
                    seen_ids.add(x)
                new_ids_all.extend(new_ids)

            if new_bps:
                append_lines(SEEN_BIOPROJECTS_PATH, sorted(new_bps))
                seen_bioprojects |= new_bps
                new_bps_all |= new_bps

            print(f"[{day.isoformat()}] +{len(added)} BioProjects, +{len(new_ids)} ids")

    else:
        # Daily: last N days up to today UTC (inclusive)
        end = dt.datetime.now(dt.timezone.utc).date()
        start = end - dt.timedelta(days=max(int(args.days), 1))
        cat_path = catalog_path_for_year(end.year)

        for day in daterange(start, end):
            try:
                added, new_ids, new_bps = run_day(
                    args.query, day, args.max_per_day, seen_ids, seen_bioprojects, bool(args.require_location)
                )
            except Exception as e:
                print(f"[{day.isoformat()}] ERROR: {e}")
                continue

            if added:
                append_jsonl(cat_path, added)
                added_all.extend(added)

            if new_ids:
                append_lines(SEEN_IDS_PATH, new_ids)
                for x in new_ids:
                    seen_ids.add(x)
                new_ids_all.extend(new_ids)

            if new_bps:
                append_lines(SEEN_BIOPROJECTS_PATH, sorted(new_bps))
                seen_bioprojects |= new_bps
                new_bps_all |= new_bps

            print(f"[{day.isoformat()}] +{len(added)} BioProjects, +{len(new_ids)} ids")

    # Public artifacts
    write_latest_json(added_all, args.query, bool(args.require_location))
    write_latest_csv(added_all)

    idx = build_archive_index()
    write_archive_html(idx)

    db_index = write_database_exports()
    write_dashboard_stub()

    print(
        f"UrbanScope done. added_bioprojects={len(added_all)} "
        f"new_bioprojects={len(new_bps_all)} new_ids={len(new_ids_all)} db_total={db_index.get('total_records')}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
