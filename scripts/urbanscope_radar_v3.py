#!/usr/bin/env python3
"""
UrbanScope — Urban Metagenomics & Metatranscriptomics Dataset Radar (GitHub-only)

What it does
------------
- Backfill year-by-year + daily incremental updates (GitHub Actions friendly)
- Finds urban metagenomics/metatranscriptomics studies and their datasets using:
    (A) SRA direct search (db=sra)
    (B) PubMed search + elink PubMed -> SRA
- Keeps a GitHub-native, continually-updated “database”:
    data/catalog_<YEAR>.jsonl        (append-only yearly catalogs)
    data/seen_ids.txt               (dedupe index for pmid:* and sra:*)
    data/seen_projects.txt          (dedupe index for bioproject PRJ*; enforces 1 SRA per project)
- Enrichment:
    - optional city/country extraction (best-effort)
    - study-type classification: air / wastewater / surface / other
- Publishes public artifacts for GitHub Pages:
    docs/index.html                 (dashboard)
    docs/db/index.json              (database index + counts)
    docs/db/year/<YEAR>.json        (year slices for fast UI)
    docs/db/records.ndjson          (full database, streamable)
    docs/db/records.csv             (full database CSV)
    docs/archive/index.html         (archive by year)
    docs/archive/index.json
    docs/archive/csv/catalog_<YEAR>.csv
    docs/archive/csv/latest_added.csv

Key behavior requested
----------------------
- City/country filter is OPTIONAL:
    --require-location  (only keep records with BOTH city and country)
- Ensure only ONE SRA is listed per BioProject (PRJNA/PRJEB/PRJDB):
    - persisted via data/seen_projects.txt
    - if BioProject cannot be determined, fallback to a per-uid pseudo-project key UID:<uid>
      (still enforces “only one” for that unknown project)

Usage
-----
  # Backfill one year
  python scripts/urbanscope_radar.py backfill-year --year 2018 --max-per-day 200

  # Daily incremental update (last N days)
  python scripts/urbanscope_radar.py daily --days 1 --max-per-day 200

  # Strict location requirement (optional)
  python scripts/urbanscope_radar.py daily --days 1 --require-location

Notes
-----
- 100% stdlib
- NCBI E-utilities sometimes disconnects: robust retry/backoff is built in
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

SEEN_PATH = os.path.join(DATA_DIR, "seen_ids.txt")
SEEN_PROJECTS_PATH = os.path.join(DATA_DIR, "seen_projects.txt")

LATEST_JSON_PATH = os.path.join(DOCS_DIR, "latest.json")
ARCHIVE_INDEX_JSON = os.path.join(ARCHIVE_DIR, "index.json")
ARCHIVE_INDEX_HTML = os.path.join(ARCHIVE_DIR, "index.html")
DOCS_INDEX_HTML = os.path.join(DOCS_DIR, "index.html")
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
                    country = country or ("United States" if c in {"USA", "United States"} else c)
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
        # Occasionally upstream returns HTML error payloads; surface URL for debugging.
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


def esummary_sra(uids: List[str]) -> List[Dict[str, Any]]:
    """
    Slightly richer SRA summary:
    - Title
    - Best-effort BioProject accession (PRJ*)
    Still lightweight enough for GitHub.
    """
    if not uids:
        return []
    params = {"db": "sra", "id": ",".join(uids), "retmode": "xml"}
    url = EUTILS + "esummary.fcgi?" + urllib.parse.urlencode(params)
    root = parse_xml(http_get(url), url)

    items: List[Dict[str, Any]] = []
    for docsum in root.findall(".//DocSum"):
        uid = (docsum.findtext("Id") or "").strip()
        rec: Dict[str, Any] = {"sra_uid": uid, "title": "", "bioproject": ""}

        for it in docsum.findall("Item"):
            name = (it.attrib.get("Name", "") or "").strip()
            text = (it.text or "").strip()

            if name == "Title":
                rec["title"] = text

            # NCBI can be inconsistent; try multiple name variants
            if name.lower() in {"bioproject", "bioprojectaccn", "bioprojectaccession", "bioprojectacc"}:
                m = BIOPROJECT_RE.search(text)
                rec["bioproject"] = (m.group(0).upper() if m else text)

        # Fallback: mine PRJ* from the title if not found
        if not rec["bioproject"] and rec["title"]:
            m = BIOPROJECT_RE.search(rec["title"])
            if m:
                rec["bioproject"] = m.group(0).upper()

        rec["sra_link"] = f"https://www.ncbi.nlm.nih.gov/sra/?term={uid}" if uid else ""
        items.append(rec)

    return items


# -------------------------
# Dedupe + storage helpers
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


def daterange(start: dt.date, end: dt.date) -> Iterable[dt.date]:
    cur = start
    while cur <= end:
        yield cur
        cur += dt.timedelta(days=1)


# -------------------------
# Core acquisition: enforce 1 SRA per BioProject
# -------------------------
def project_key_from_summary(s: Dict[str, Any]) -> str:
    """
    Returns a stable key for "one SRA per project".
    Prefer BioProject accession; fallback to UID pseudo-key.
    """
    bp = (s.get("bioproject", "") or "").strip().upper()
    uid = (s.get("sra_uid", "") or "").strip()
    if bp:
        return bp
    # fallback: treat UID as its own project bucket
    return f"UID:{uid}" if uid else ""


def run_day(
    query: str,
    day: dt.date,
    retmax: int,
    seen_ids: Set[str],
    require_location: bool,
    seen_projects: Set[str],
) -> Tuple[List[Dict[str, Any]], List[str], Set[int], Set[str]]:
    """
    Returns:
      new_records, new_seen_ids, touched_years, new_projects_added
    """
    mindate = day.strftime("%Y/%m/%d")
    maxdate = day.strftime("%Y/%m/%d")

    new_records: List[Dict[str, Any]] = []
    new_seen: List[str] = []
    touched_years: Set[int] = {day.year}
    new_projects: Set[str] = set()

    # Polite throttling between calls (helps avoid disconnects)
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
        try:
            summaries = esummary_sra(sra_uids_new[:300])
        except Exception as e:
            print(f"[{day.isoformat()}] SRA esummary ERROR: {e}")
            summaries = []

        for s in summaries:
            uid = (s.get("sra_uid") or "").strip()
            title = (s.get("title") or "").strip()
            proj_key = project_key_from_summary(s)

            if not uid or not proj_key:
                continue

            # Enforce one SRA per project
            if proj_key in seen_projects or proj_key in new_projects:
                continue

            stype = classify_study_type(title)
            city, country = extract_city_country(title)

            if require_location and (not city or not country):
                continue

            new_projects.add(proj_key)

            new_records.append(
                {
                    "source": "sra_search",
                    "day_utc": day.isoformat(),
                    "query": query,
                    "study_type": stype,
                    "city": city,
                    "country": country,
                    "title": title,
                    "sra_uid": uid,
                    "bioproject": (s.get("bioproject") or "").strip().upper() or None,
                    "sra_link": s.get("sra_link", ""),
                    "ingested_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                }
            )
            new_seen.append(f"sra:{uid}")

    nap()

    # -----------------
    # B) PubMed -> SRA (choose 1 representative SRA per project)
    # -----------------
    try:
        pmids = esearch("pubmed", query, mindate, maxdate, retmax=retmax, sort="pub+date")
    except Exception as e:
        print(f"[{day.isoformat()}] PubMed esearch ERROR: {e}")
        pmids = []
    nap()

    pmids_new = [p for p in pmids if f"pmid:{p}" not in seen_ids]
    if pmids_new:
        try:
            pmid_to_sra = elink_pubmed_to_sra(pmids_new[:200])
        except Exception as e:
            print(f"[{day.isoformat()}] PubMed elink ERROR: {e}")
            pmid_to_sra = {}

        # Flatten uids and summarize once so we can pick 1 per project reliably
        all_uids = sorted({uid for uids in pmid_to_sra.values() for uid in uids})
        uid_summary: Dict[str, Dict[str, Any]] = {}
        if all_uids:
            try:
                sums = esummary_sra(all_uids[:400])
                uid_summary = {str(r.get("sra_uid")): r for r in sums if r.get("sra_uid")}
            except Exception as e:
                print(f"[{day.isoformat()}] SRA esummary (from PubMed links) ERROR: {e}")
                uid_summary = {}

        nap()

        for pmid, uids in pmid_to_sra.items():
            if not uids:
                continue

            chosen_uid: Optional[str] = None
            chosen_summary: Dict[str, Any] = {}

            # Select first UID whose project has not been seen
            for uid in uids:
                s = uid_summary.get(uid, {"sra_uid": uid, "title": "", "bioproject": ""})
                proj_key = project_key_from_summary(s)
                if not proj_key:
                    continue
                if proj_key in seen_projects or proj_key in new_projects:
                    continue
                chosen_uid = uid
                chosen_summary = s
                new_projects.add(proj_key)
                break

            if not chosen_uid:
                continue

            title = (chosen_summary.get("title") or "").strip()
            stype = classify_study_type(title or query)
            city, country = extract_city_country(title) if title else (None, None)

            if require_location and (not city or not country):
                # If strict mode is on and we couldn't extract location, skip.
                continue

            new_records.append(
                {
                    "source": "pubmed_elink_sra",
                    "day_utc": day.isoformat(),
                    "query": query,
                    "study_type": stype,
                    "city": city,
                    "country": country,
                    "title": title or None,
                    "pmid": pmid,
                    "paper": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                    "sra_uid": chosen_uid,
                    "bioproject": (chosen_summary.get("bioproject") or "").strip().upper() or None,
                    "sra_link": f"https://www.ncbi.nlm.nih.gov/sra/?term={chosen_uid}",
                    "ingested_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                }
            )
            new_seen.append(f"pmid:{pmid}")
            new_seen.append(f"sra:{chosen_uid}")

    # De-dupe IDs within-run and against seen
    uniq_ids: List[str] = []
    local = set()
    for x in new_seen:
        if x not in seen_ids and x not in local:
            uniq_ids.append(x)
            local.add(x)

    # Filter any direct SRA records that became "seen" mid-run
    filtered_records: List[Dict[str, Any]] = []
    for r in new_records:
        if r.get("source") == "sra_search":
            uid = r.get("sra_uid", "")
            if uid and f"sra:{uid}" in seen_ids:
                continue
        filtered_records.append(r)

    return filtered_records, uniq_ids, touched_years, new_projects


# -------------------------
# Database exports + archive + dashboard
# -------------------------
CSV_FIELDS = [
    "source",
    "day_utc",
    "study_type",
    "city",
    "country",
    "title",
    "pmid",
    "paper",
    "sra_uid",
    "bioproject",
    "sra_link",
    "ingested_utc",
]

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


def write_latest_json(records: List[Dict[str, Any]], query: str, require_location: bool) -> None:
    payload = {
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "query": query,
        "require_location": bool(require_location),
        "added_records": len(records),
        "records": records[-300:],
    }
    with open(LATEST_JSON_PATH, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)


def write_latest_csv(records: List[Dict[str, Any]]) -> None:
    with open(LATEST_CSV_PATH, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        w.writeheader()
        for r in records[-5000:]:
            w.writerow({k: r.get(k, "") for k in CSV_FIELDS})


def write_year_csv(year: int) -> Optional[str]:
    src = catalog_path_for_year(year)
    if not os.path.exists(src):
        return None
    out = os.path.join(ARCHIVE_CSV_DIR, f"catalog_{year}.csv")
    with open(out, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        w.writeheader()
        for r in iter_jsonl(src):
            w.writerow({k: r.get(k, "") for k in CSV_FIELDS})
    return out


def build_archive_index() -> Dict[str, Any]:
    years = []
    for y in list_catalog_years():
        path = catalog_path_for_year(y)
        count = 0
        with open(path, "r", encoding="utf-8") as f:
            for _ in f:
                count += 1
        years.append({"year": y, "records": count, "csv": f"csv/catalog_{y}.csv"})
    payload = {
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "years": years,
    }
    with open(ARCHIVE_INDEX_JSON, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)
    return payload


def write_archive_html(index_payload: Dict[str, Any]) -> None:
    years = index_payload.get("years", [])
    rows = []
    for y in years:
        rows.append(
            f"<tr><td><b>{y['year']}</b></td><td>{y['records']}</td>"
            f"<td><a href='{y['csv']}'>catalog_{y['year']}.csv</a></td></tr>"
        )
    html = f"""<!doctype html>
<html><head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1"/>
<title>UrbanScope Archive</title>
<style>
body{{font-family:system-ui,Segoe UI,Roboto,Arial;max-width:980px;margin:40px auto;padding:0 16px;line-height:1.55}}
.muted{{color:#666}}
table{{border-collapse:collapse;width:100%;margin-top:14px}}
th,td{{border-bottom:1px solid #eee;padding:10px;text-align:left}}
code{{background:#f2f2f2;padding:2px 6px;border-radius:6px}}
</style>
</head><body>
<h1>UrbanScope Archive</h1>
<p class="muted">Yearly catalogs exported as CSV for collaborators.</p>
<p class="muted">Index updated: <code>{index_payload.get("generated_at_utc","")}</code></p>
<p><a href="../index.html">← Back to dashboard</a> · <a href="../db/index.json">db/index.json</a></p>
<table>
<thead><tr><th>Year</th><th>Records</th><th>CSV</th></tr></thead>
<tbody>
{''.join(rows) if rows else '<tr><td colspan="3" class="muted">No archive yet.</td></tr>'}
</tbody>
</table>
</body></html>
"""
    with open(ARCHIVE_INDEX_HTML, "w", encoding="utf-8") as f:
        f.write(html)


def write_database_exports() -> Dict[str, Any]:
    years = list_catalog_years()
    total = 0
    by_year: Dict[int, int] = {}
    by_type: Dict[str, int] = {}
    by_country: Dict[str, int] = {}

    # Full NDJSON
    with open(DB_ALL_NDJSON, "w", encoding="utf-8") as nd:
        for y in years:
            for r in iter_jsonl(catalog_path_for_year(y)):
                nd.write(json.dumps(r, ensure_ascii=False) + "\n")
                total += 1
                by_year[y] = by_year.get(y, 0) + 1
                st = r.get("study_type") or "other"
                by_type[st] = by_type.get(st, 0) + 1
                c = r.get("country") or "Unknown"
                by_country[c] = by_country.get(c, 0) + 1

    # Full CSV
    with open(DB_ALL_CSV, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        w.writeheader()
        for y in years:
            for r in iter_jsonl(catalog_path_for_year(y)):
                w.writerow({k: r.get(k, "") for k in CSV_FIELDS})

    # Per-year JSON slices (for fast UI)
    for y in years:
        out = os.path.join(DB_YEAR_DIR, f"{y}.json")
        items = list(iter_jsonl(catalog_path_for_year(y)))
        with open(out, "w", encoding="utf-8") as f:
            json.dump({"year": y, "records": items}, f, ensure_ascii=False)

    index = {
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "years": years,
        "total_records": total,
        "counts_by_year": {str(k): v for k, v in sorted(by_year.items())},
        "counts_by_type": {k: v for k, v in sorted(by_type.items(), key=lambda x: x[0])},
        "top_countries": sorted(by_country.items(), key=lambda x: x[1], reverse=True)[:25],
        "downloads": {
            "csv": "db/records.csv",
            "ndjson": "db/records.ndjson",
            "archive": "archive/index.html",
        },
    }
    with open(DB_INDEX_JSON, "w", encoding="utf-8") as f:
        json.dump(index, f, indent=2, ensure_ascii=False)

    return index


def write_dashboard_html() -> None:
    """docs/index.html – DB browser (year slices) + latest-run view + filtering."""
    html = """<!doctype html>
<html><head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1"/>
<title>UrbanScope</title>
<style>
body{font-family:system-ui,Segoe UI,Roboto,Arial;max-width:1100px;margin:40px auto;padding:0 16px;line-height:1.55}
.muted{color:#666}
a{text-decoration:none} a:hover{text-decoration:underline}
.tag{display:inline-block;background:#f3f3f3;border-radius:999px;padding:3px 10px;font-size:.85em;color:#444;margin-left:8px}
.controls{display:flex;gap:10px;flex-wrap:wrap;margin:16px 0}
input,select{padding:10px;border:1px solid #ddd;border-radius:10px}
.card{border:1px solid #eee;border-radius:14px;padding:14px;margin:10px 0}
code{background:#f2f2f2;padding:2px 6px;border-radius:6px}
.grid{display:grid;grid-template-columns:1fr 1fr;gap:14px}
@media(max-width:900px){.grid{grid-template-columns:1fr}}
</style>
</head><body>
<h1>UrbanScope <span class="tag">urban metagenomics • metatranscriptomics</span></h1>

<div class="grid">
  <div>
    <h2 style="margin:0 0 6px 0;">Database</h2>
    <p class="muted" id="dbmeta">Loading db/index.json…</p>
    <p>
      <a href="db/records.csv">Download full CSV</a> ·
      <a href="db/records.ndjson">Download full NDJSON</a> ·
      <a href="archive/index.html">Browse archive by year</a>
    </p>
  </div>
  <div>
    <h2 style="margin:0 0 6px 0;">Latest run</h2>
    <p class="muted" id="runmeta">Loading latest.json…</p>
    <p><a href="archive/csv/latest_added.csv">Download latest_added.csv</a></p>
  </div>
</div>

<div class="controls">
  <select id="year"></select>
  <select id="stype">
    <option value="">All study types</option>
    <option value="wastewater">wastewater</option>
    <option value="air">air</option>
    <option value="surface">surface</option>
    <option value="other">other</option>
  </select>
  <input id="city" placeholder="City filter" />
  <input id="country" placeholder="Country filter" />
  <input id="q" placeholder="Search title / links…" style="flex:1;min-width:220px"/>
</div>

<div id="list"></div>

<script>
function norm(x){ return (x||"").toString().toLowerCase(); }

let DB = {years:[]};
let YEAR_DATA = {records:[]};

async function loadDB(){
  const res = await fetch("db/index.json", {cache:"no-store"});
  DB = await res.json();
  const years = (DB.years || []).map(String);
  const ysel = document.getElementById("year");
  ysel.innerHTML = years.map(y => `<option value="${y}">${y}</option>`).join("");
  const latestYear = years.length ? years[years.length-1] : "";
  if (latestYear) ysel.value = latestYear;

  document.getElementById("dbmeta").textContent =
    `Total records: ${DB.total_records} • Years: ${(DB.years||[]).join(", ")} • Updated (UTC): ${DB.generated_at_utc}`;

  if (latestYear) await loadYear(latestYear);
}

async function loadRun(){
  const res = await fetch("latest.json", {cache:"no-store"});
  const run = await res.json();
  document.getElementById("runmeta").textContent =
    `Added records: ${run.added_records} • require_location: ${run.require_location} • Updated (UTC): ${run.generated_at_utc}`;
}

async function loadYear(y){
  const res = await fetch(`db/year/${y}.json`, {cache:"no-store"});
  YEAR_DATA = await res.json();
  render();
}

function render(){
  const records = YEAR_DATA.records || [];
  const q = norm(document.getElementById("q").value);
  const st = norm(document.getElementById("stype").value);
  const city = norm(document.getElementById("city").value);
  const country = norm(document.getElementById("country").value);

  const filtered = records.filter(r => {
    const blob = norm((r.title||"") + " " + (r.paper||"") + " " + (r.sra_link||"") + " " + (r.sra_uid||""));
    const okQ = !q || blob.includes(q);
    const okS = !st || norm(r.study_type) === st;
    const okC = !city || norm(r.city).includes(city);
    const okK = !country || norm(r.country).includes(country);
    return okQ && okS && okC && okK;
  });

  const list = document.getElementById("list");
  if (!filtered.length) {
    list.innerHTML = "<p class='muted'>No records match filters for this year.</p>";
    return;
  }

  list.innerHTML = filtered.slice().reverse().map(r => {
    const src = r.source || "";
    const st = r.study_type ? `<span class="tag">${r.study_type}</span>` : "";
    const loc = [r.city, r.country].filter(Boolean).join(", ");
    const locHtml = loc ? `<div class="muted" style="margin-top:6px">📍 ${loc}</div>` : "";
    const title = r.title ? `<div style="margin-top:8px">${r.title}</div>` : "";
    const paper = r.paper ? `<div><a href="${r.paper}" target="_blank" rel="noreferrer"><b>Paper</b></a></div>` : "";
    const dataset = r.sra_link ? `<div><a href="${r.sra_link}" target="_blank" rel="noreferrer"><b>Dataset</b></a></div>` : "";
    const uid = r.sra_uid ? `<div class="muted">SRA UID: <code>${r.sra_uid}</code></div>` : "";
    const bp = r.bioproject ? `<div class="muted">BioProject: <code>${r.bioproject}</code></div>` : "";
    const pmid = r.pmid ? `<div class="muted">PMID: <code>${r.pmid}</code></div>` : "";
    return `<div class="card">
      <div class="muted">${src} ${st}</div>
      ${paper}
      ${dataset}
      ${uid}
      ${bp}
      ${pmid}
      ${locHtml}
      ${title}
    </div>`;
  }).join("");
}

document.getElementById("year").addEventListener("change", async (e) => loadYear(e.target.value));
["q","stype","city","country"].forEach(id => {
  document.getElementById(id).addEventListener(id==="stype" ? "change" : "input", render);
});

loadDB();
loadRun();
</script>
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

    seen_ids = load_lineset(SEEN_PATH)
    seen_projects = load_lineset(SEEN_PROJECTS_PATH)

    added_records_all: List[Dict[str, Any]] = []
    touched_years: Set[int] = set()
    added_ids_all: List[str] = []

    if args.cmd == "backfill-year":
        year = int(args.year)
        start = dt.date(year, 1, 1)
        end = dt.date(year, 12, 31)
        cat_path = catalog_path_for_year(year)

        for day in daterange(start, end):
            try:
                recs, new_ids, tys, new_projects = run_day(
                    args.query, day, args.max_per_day, seen_ids, args.require_location, seen_projects
                )
            except Exception as e:
                print(f"[{day.isoformat()}] ERROR: {e}")
                continue

            if recs:
                append_jsonl(cat_path, recs)
                added_records_all.extend(recs)

            if new_ids:
                append_lines(SEEN_PATH, new_ids)
                for x in new_ids:
                    seen_ids.add(x)
                added_ids_all.extend(new_ids)

            if new_projects:
                append_lines(SEEN_PROJECTS_PATH, sorted(new_projects))
                seen_projects |= new_projects

            touched_years |= tys
            print(f"[{day.isoformat()}] +{len(recs)} records, +{len(new_ids)} ids, +{len(new_projects)} projects")

    else:
        # Daily: last N days up to today UTC
        end = dt.datetime.now(dt.timezone.utc).date()
        start = end - dt.timedelta(days=max(int(args.days), 1))
        cat_path = catalog_path_for_year(end.year)

        for day in daterange(start, end):
            try:
                recs, new_ids, tys, new_projects = run_day(
                    args.query, day, args.max_per_day, seen_ids, args.require_location, seen_projects
                )
            except Exception as e:
                print(f"[{day.isoformat()}] ERROR: {e}")
                continue

            if recs:
                append_jsonl(cat_path, recs)
                added_records_all.extend(recs)

            if new_ids:
                append_lines(SEEN_PATH, new_ids)
                for x in new_ids:
                    seen_ids.add(x)
                added_ids_all.extend(new_ids)

            if new_projects:
                append_lines(SEEN_PROJECTS_PATH, sorted(new_projects))
                seen_projects |= new_projects

            touched_years |= tys
            print(f"[{day.isoformat()}] +{len(recs)} records, +{len(new_ids)} ids, +{len(new_projects)} projects")

    # Views / exports
    write_latest_json(added_records_all, args.query, bool(args.require_location))
    write_latest_csv(added_records_all)

    # Update yearly CSVs for touched years (fast)
    for y in sorted(touched_years):
        out_csv = write_year_csv(y)
        if out_csv:
            print(f"Wrote {out_csv}")

    # Archive index + html
    idx = build_archive_index()
    write_archive_html(idx)

    # Database exports + dashboard
    db_index = write_database_exports()
    write_dashboard_html()

    print(
        f"UrbanScope done. added_records={len(added_records_all)} "
        f"added_ids={len(added_ids_all)} touched_years={sorted(touched_years) if touched_years else []} "
        f"db_total={db_index.get('total_records')}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
