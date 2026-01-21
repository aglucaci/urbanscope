#!/usr/bin/env python3
"""
UrbanScope — Urban Metagenomics & Metatranscriptomics Dataset Radar (GitHub-only)

What it does
------------
- Historical backfill (year-by-year) + daily incremental updates (GitHub Actions-friendly)
- Finds URBAN metagenomics / metatranscriptomics studies and links them to datasets:
  (A) SRA direct search (db=sra)
  (B) PubMed search + elink PubMed -> SRA (db=sra)
- Permanent dedupe across runs via data/seen_ids.txt
- Append-only yearly catalogs: data/catalog_<YEAR>.jsonl
- Adds lightweight enrichment:
  - City / country extraction (heuristics from titles)
  - Study-type classification: air / wastewater / surface / other
- Publishes GitHub Pages-ready artifacts:
  - docs/latest.json (new records added this run)
  - docs/index.html (dashboard)
  - docs/archive/index.json + docs/archive/index.html (browse years)
  - docs/archive/csv/catalog_<YEAR>.csv (collaborator-friendly exports per year)
  - docs/archive/csv/latest_added.csv (just-added records as CSV)

Run modes
---------
  # Backfill one year (recommended for GitHub Actions)
  python scripts/urbanscope_radar.py backfill-year --year 2016

  # Daily (last N days up to today UTC, writes into current year's catalog)
  python scripts/urbanscope_radar.py daily --days 1

Notes
-----
- 100% stdlib (no dependencies)
- City/country extraction is heuristic (title-only); it will be "best effort".
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
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple
import random
import socket
import http.client
import urllib.error

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

DEFAULT_QUERY = (
    '(urban OR city OR cities OR metropolis OR megacity OR subway OR transit OR "public transit" '
    'OR "built environment" OR "indoor microbiome" OR wastewater OR sewage OR stormwater) '
    'AND (metagenomics OR metatranscriptomics OR "shotgun metagenomic" OR "RNA-seq" OR "total RNA")'
)

# -------------------------
# Paths / folders (GitHub-only)
# -------------------------
DATA_DIR = "data"
DOCS_DIR = "docs"
ARCHIVE_DIR = os.path.join(DOCS_DIR, "archive")
ARCHIVE_CSV_DIR = os.path.join(ARCHIVE_DIR, "csv")

SEEN_PATH = os.path.join(DATA_DIR, "seen_ids.txt")
LATEST_JSON_PATH = os.path.join(DOCS_DIR, "latest.json")
LATEST_CSV_PATH = os.path.join(ARCHIVE_CSV_DIR, "latest_added.csv")
ARCHIVE_INDEX_JSON = os.path.join(ARCHIVE_DIR, "index.json")
ARCHIVE_INDEX_HTML = os.path.join(ARCHIVE_DIR, "index.html")
DOCS_INDEX_HTML = os.path.join(DOCS_DIR, "index.html")
DB_DIR = os.path.join(DOCS_DIR, "db")
DB_YEAR_DIR = os.path.join(DB_DIR, "year")
DB_INDEX_JSON = os.path.join(DB_DIR, "index.json")
DB_ALL_NDJSON = os.path.join(DB_DIR, "records.ndjson")
DB_ALL_CSV = os.path.join(DB_DIR, "records.csv")

os.makedirs(DB_DIR, exist_ok=True)
os.makedirs(DB_YEAR_DIR, exist_ok=True)
# -------------------------
# Country / city heuristics
# -------------------------
# A compact but useful country list (common + ISO-ish). You can extend later.
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
    key=lambda s: (-len(s), s.lower())
)

# A small city alias table (high-signal urban omics locales). Extend as you like.
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

# Study-type keyword rules
STUDY_TYPE_RULES = [
    ("wastewater", [r"\bwastewater\b", r"\bsewage\b", r"\bstormwater\b", r"\binfluent\b", r"\beffluent\b", r"\bWWTP\b"]),
    ("air", [r"\bair\b", r"\bairborne\b", r"\baerosol\b", r"\baerosols\b", r"\bHVAC\b", r"\bventilation\b", r"\bdust\b"]),
    ("surface", [r"\bsurface\b", r"\bfomite\b", r"\btouch\b", r"\bhandrail\b", r"\bswab\b", r"\bdoorknob\b", r"\bbench\b"]),
]


# -------------------------
# HTTP / NCBI helpers
# -------------------------
def http_get(url: str, timeout: int = 30, retries: int = 6, base_sleep: float = 0.8) -> bytes:
    """
    Robust GET with retries/backoff for flaky upstream (e.g., NCBI E-utilities).

    Retries on:
      - RemoteDisconnected
      - timeout / socket errors
      - HTTP 429/500/502/503/504
      - transient URLError
    """
    last_err = None

    for attempt in range(retries):
        try:
            # NCBI asks for identifying UA; keep it stable.
            req = urllib.request.Request(
                url,
                headers={
                    "User-Agent": "urbanscope/1.0 (GitHub Actions; contact: none)",
                    "Accept": "application/xml,text/xml,application/json;q=0.9,*/*;q=0.8",
                },
            )
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = resp.read()
                # Small guard: if upstream returns empty, treat as retryable
                if not data:
                    raise urllib.error.URLError("Empty response body")
                return data

        except (http.client.RemoteDisconnected, TimeoutError, socket.timeout, ConnectionResetError) as e:
            last_err = e

        except urllib.error.HTTPError as e:
            # Retry on transient server/rate-limit errors
            if e.code in (429, 500, 502, 503, 504):
                last_err = e
            else:
                raise

        except urllib.error.URLError as e:
            # Often transient network issues; retry
            last_err = e

        # Exponential backoff + jitter
        sleep_s = base_sleep * (2 ** attempt) + random.uniform(0, 0.25)
        time.sleep(sleep_s)

    raise RuntimeError(f"HTTP GET failed after {retries} retries. Last error: {last_err}")


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
    """Lightweight SRA summary (Title only) to avoid huge repo artifacts."""
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
            if name == "Title":
                rec["title"] = text
        rec["sra_link"] = f"https://www.ncbi.nlm.nih.gov/sra/?term={uid}" if uid else ""
        items.append(rec)
    return items

def list_catalog_years() -> List[int]:
    years = []
    for fn in os.listdir(DATA_DIR):
        m = re.match(r"catalog_(\d{4})\.jsonl$", fn)
        if m:
            years.append(int(m.group(1)))
    return sorted(years)

def read_all_records(years: List[int]) -> Iterable[Dict[str, Any]]:
    for y in years:
        yield from iter_jsonl(catalog_path_for_year(y))

def write_database_exports() -> Dict[str, Any]:
    years = list_catalog_years()

    total = 0
    by_year = {}
    by_type = {}
    by_country = {}

    # Write NDJSON (streamable full DB)
    with open(DB_ALL_NDJSON, "w", encoding="utf-8") as nd:
        for r in read_all_records(years):
            nd.write(json.dumps(r, ensure_ascii=False) + "\n")
            total += 1
            yr = int((r.get("day_utc") or "0000")[:4] or 0)
            by_year[yr] = by_year.get(yr, 0) + 1
            st = r.get("study_type") or "other"
            by_type[st] = by_type.get(st, 0) + 1
            c = r.get("country") or "Unknown"
            by_country[c] = by_country.get(c, 0) + 1

    # Write CSV (full DB)
    with open(DB_ALL_CSV, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        w.writeheader()
        for r in read_all_records(years):
            row = {k: r.get(k, "") for k in CSV_FIELDS}
            if isinstance(row.get("sra_uids"), list):
                row["sra_uids"] = " ".join(row["sra_uids"])
            w.writerow(row)

    # Optional: per-year JSON slices for faster front-end loads
    for y in years:
        out = os.path.join(DB_YEAR_DIR, f"{y}.json")
        items = list(iter_jsonl(catalog_path_for_year(y)))
        with open(out, "w", encoding="utf-8") as f:
            json.dump({"year": y, "records": items}, f, ensure_ascii=False)

    # Index metadata
    index = {
        "generated_at_utc": dt.datetime.utcnow().isoformat(timespec="seconds"),
        "years": years,
        "total_records": total,
        "counts_by_year": dict(sorted(by_year.items())),
        "counts_by_type": dict(sorted(by_type.items(), key=lambda x: x[0])),
        # top 25 countries by count
        "top_countries": sorted(by_country.items(), key=lambda x: x[1], reverse=True)[:25],
        "downloads": {
            "ndjson": "db/records.ndjson",
            "csv": "db/records.csv",
            "archive": "archive/index.html",
        },
    }

    with open(DB_INDEX_JSON, "w", encoding="utf-8") as f:
        json.dump(index, f, indent=2, ensure_ascii=False)

    return index

# -------------------------
# Dedupe + storage
# -------------------------
def ensure_dirs() -> None:
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(DOCS_DIR, exist_ok=True)
    os.makedirs(ARCHIVE_DIR, exist_ok=True)
    os.makedirs(ARCHIVE_CSV_DIR, exist_ok=True)


def load_seen(path: str) -> Set[str]:
    if not os.path.exists(path):
        return set()
    with open(path, "r", encoding="utf-8") as f:
        return set(line.strip() for line in f if line.strip())


def append_seen(path: str, ids: Iterable[str]) -> None:
    ids = list(ids)
    if not ids:
        return
    with open(path, "a", encoding="utf-8") as f:
        for i in ids:
            f.write(i + "\n")


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
# Enrichment: study type, city, country
# -------------------------
def classify_study_type(text: str) -> str:
    t = (text or "").lower()
    for label, patterns in STUDY_TYPE_RULES:
        for pat in patterns:
            if re.search(pat, t, flags=re.IGNORECASE):
                return label
    return "other"


def extract_city_country(text: str) -> Tuple[Optional[str], Optional[str]]:
    """Best-effort city/country extraction from title-like text."""
    if not text:
        return None, None

    # Country: scan for known names first (longest first)
    country = None
    for c in COUNTRIES:
        if re.search(rf"\b{re.escape(c)}\b", text, flags=re.IGNORECASE):
            # Normalize a few common aliases
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

    # City: match alias patterns
    city = None
    for canonical, patterns in CITY_ALIASES.items():
        for pat in patterns:
            if re.search(pat, text, flags=re.IGNORECASE):
                city = canonical
                break
        if city:
            break

    # Heuristic: "City, Country" pattern
    if not city:
        m = re.search(r"\b([A-Z][A-Za-z.\- ]{2,40}),\s*([A-Z][A-Za-z.\- ]{2,40})\b", text)
        if m:
            cand_city = m.group(1).strip()
            cand_country = m.group(2).strip()
            # If the second token looks like a known country, prefer it
            for c in COUNTRIES:
                if cand_country.lower() == c.lower() or re.search(rf"\b{re.escape(c)}\b", cand_country, flags=re.I):
                    country = country or ( "United States" if c in {"USA", "United States"} else c )
                    break
            # Accept city if it looks plausible (not too long)
            if 2 <= len(cand_city) <= 40:
                city = cand_city

    return city, country


# -------------------------
# Core acquisition logic
# -------------------------
def run_day(query: str, day: dt.date, retmax: int, seen: Set[str]) -> Tuple[List[Dict[str, Any]], List[str], Set[int]]:
    """
    For a given day window, collect:
      - SRA hits directly (db=sra)
      - PubMed hits linked to SRA via elink
    Returns:
      new_records, new_seen_ids, touched_years (years that should get CSV regenerated)
    """
    mindate = day.strftime("%Y/%m/%d")
    maxdate = day.strftime("%Y/%m/%d")

    new_records: List[Dict[str, Any]] = []
    new_seen: List[str] = []
    touched_years: Set[int] = {day.year}

    # A) SRA direct
    try:
        sra_uids = esearch("sra", query, mindate, maxdate, retmax=retmax, sort="Date")
    except Exception:
        sra_uids = []
    time.sleep(0.6)

    sra_uids_new = [u for u in sra_uids if f"sra:{u}" not in seen]
    if sra_uids_new:
        summaries = esummary_sra(sra_uids_new[:200])
        for s in summaries:
            uid = s.get("sra_uid", "")
            title = s.get("title", "") or ""
            stype = classify_study_type(title)
            city, country = extract_city_country(title)
            # Enforce location requirement
            if not city or not country:
                continue 
            if uid:
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
                        "sra_link": s.get("sra_link", ""),
                        "ingested_utc": dt.datetime.utcnow().isoformat(timespec="seconds"),
                    }
                )
                new_seen.append(f"sra:{uid}")
    time.sleep(0.34)

    # B) PubMed -> SRA
    try:
        pmids = esearch("pubmed", query, mindate, maxdate, retmax=retmax, sort="pub+date")
    except Exception:
        pmids = []
    time.sleep(0.34)

    pmids_new = [p for p in pmids if f"pmid:{p}" not in seen]
    if pmids_new:
        pmid_to_sra = elink_pubmed_to_sra(pmids_new[:200])
        for pmid, uids in pmid_to_sra.items():
            # We don't have the PubMed title in this minimal pipeline; classify by query/day only.
            # You can enrich later by efetch PubMed titles if you want.
            stype = classify_study_type(query)
            city, country = extract_city_country(query)
            # Enforce: must have BOTH city and country
            if not city or not country:
                continue
            new_records.append(
                {
                    "source": "pubmed_elink_sra",
                    "day_utc": day.isoformat(),
                    "query": query,
                    "study_type": stype,
                    "city": city,
                    "country": country,
                    "pmid": pmid,
                    "paper": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                    "sra_uids": uids,
                    "sra_link": f"https://www.ncbi.nlm.nih.gov/sra/?term={' OR '.join(uids)}" if uids else "",
                    "ingested_utc": dt.datetime.utcnow().isoformat(timespec="seconds"),
                }
            )
            new_seen.append(f"pmid:{pmid}")
            for uid in uids:
                new_seen.append(f"sra:{uid}")
    time.sleep(0.34)

    # De-dupe IDs within-run & against seen
    uniq_ids: List[str] = []
    local = set()
    for x in new_seen:
        if x not in seen and x not in local:
            uniq_ids.append(x)
            local.add(x)

    # Filter direct SRA records that might now be "seen" (extra safety)
    filtered_records: List[Dict[str, Any]] = []
    for r in new_records:
        if r.get("source") == "sra_search":
            uid = r.get("sra_uid", "")
            if uid and f"sra:{uid}" in seen:
                continue
        filtered_records.append(r)

    return filtered_records, uniq_ids, touched_years


# -------------------------
# CSV export + archive generation
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
    "sra_link",
    "sra_uids",
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


def write_year_csv(year: int) -> Optional[str]:
    """Generate docs/archive/csv/catalog_<YEAR>.csv from data/catalog_<YEAR>.jsonl."""
    src = catalog_path_for_year(year)
    if not os.path.exists(src):
        return None

    out_path = os.path.join(ARCHIVE_CSV_DIR, f"catalog_{year}.csv")
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        w.writeheader()
        for r in iter_jsonl(src):
            row = {k: r.get(k, "") for k in CSV_FIELDS}
            # Flatten list of sra_uids (for pubmed_elink_sra)
            if isinstance(row.get("sra_uids"), list):
                row["sra_uids"] = " ".join(row["sra_uids"])
            w.writerow(row)

    return out_path


def build_archive_index() -> Dict[str, Any]:
    """Scan data/catalog_*.jsonl and build docs/archive/index.json with counts."""
    years = []
    for fn in os.listdir(DATA_DIR):
        m = re.match(r"catalog_(\d{4})\.jsonl$", fn)
        if not m:
            continue
        year = int(m.group(1))
        path = os.path.join(DATA_DIR, fn)
        # line-count is fast and good enough
        count = 0
        with open(path, "r", encoding="utf-8") as f:
            for _ in f:
                count += 1
        years.append({"year": year, "records": count, "csv": f"csv/catalog_{year}.csv"})
    years.sort(key=lambda x: x["year"])

    payload = {
        "generated_at_utc": dt.datetime.utcnow().isoformat(timespec="seconds"),
        "years": years,
    }

    with open(ARCHIVE_INDEX_JSON, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)

    return payload


def write_archive_html(index_payload: Dict[str, Any]) -> None:
    years = index_payload.get("years", [])
    rows = []
    for y in years:
        yr = y["year"]
        rec = y["records"]
        csv_rel = y["csv"]
        rows.append(
            f"<tr><td><b>{yr}</b></td><td>{rec}</td>"
            f"<td><a href='{csv_rel}'>catalog_{yr}.csv</a></td></tr>"
        )

    html = f"""<!doctype html>
<html>
<head>
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
</head>
<body>
  <h1>UrbanScope Archive</h1>
  <p class="muted">Yearly catalogs exported as CSV for collaborators.</p>
  <p class="muted">Index updated: <code>{index_payload.get("generated_at_utc","")}</code></p>
  <p><a href="../index.html">← Back to dashboard</a></p>
  <table>
    <thead><tr><th>Year</th><th>Records</th><th>CSV</th></tr></thead>
    <tbody>
      {''.join(rows) if rows else '<tr><td colspan="3" class="muted">No archive yet.</td></tr>'}
    </tbody>
  </table>
</body>
</html>
"""
    with open(ARCHIVE_INDEX_HTML, "w", encoding="utf-8") as f:
        f.write(html)


def write_dashboard_html() -> None:
    """docs/index.html reads docs/latest.json and links to archive."""
    html = """<!doctype html>
<html>
<head>
  <meta charset="utf-8"/>
  <meta name="viewport" content="width=device-width,initial-scale=1"/>
  <title>UrbanScope</title>
  <style>
    body{font-family:system-ui,Segoe UI,Roboto,Arial;max-width:980px;margin:40px auto;padding:0 16px;line-height:1.55}
    .muted{color:#666}
    .controls{display:flex;gap:10px;flex-wrap:wrap;margin:16px 0}
    input,select{padding:10px;border:1px solid #ddd;border-radius:10px}
    .card{border:1px solid #eee;border-radius:14px;padding:14px;margin:10px 0}
    a{text-decoration:none} a:hover{text-decoration:underline}
    code{background:#f2f2f2;padding:2px 6px;border-radius:6px}
    .tag{display:inline-block;background:#f3f3f3;border-radius:999px;padding:3px 10px;font-size:.85em;margin-left:8px;color:#444}
  </style>
</head>
<body>
  <h1>UrbanScope <span class="tag">urban metagenomics • metatranscriptomics</span></h1>
  <p class="muted">Daily dataset radar (PubMed ↔ SRA) with best-effort city/country + study-type classification.</p>
  <p class="muted" id="meta">Loading…</p>
  <p>
    <a href="archive/index.html">Browse archive by year</a> ·
    <a href="archive/csv/latest_added.csv">Download latest_added.csv</a>
  </p>

  <div class="controls">
    <input id="q" placeholder="Search titles/links…" style="flex:1;min-width:220px"/>
    <select id="stype">
      <option value="">All study types</option>
      <option value="wastewater">wastewater</option>
      <option value="air">air</option>
      <option value="surface">surface</option>
      <option value="other">other</option>
    </select>
    <input id="city" placeholder="City filter" style="min-width:180px"/>
    <input id="country" placeholder="Country filter" style="min-width:180px"/>
  </div>

  <div id="list"></div>

<script>
(async () => {
  const res = await fetch("latest.json", {cache:"no-store"});
  const data = await res.json();
  const meta = document.getElementById("meta");
  meta.textContent = `Updated (UTC): ${data.generated_at_utc} • Added records this run: ${data.added_records}`;

  const records = data.records || [];
  const list = document.getElementById("list");

  const qEl = document.getElementById("q");
  const sEl = document.getElementById("stype");
  const cEl = document.getElementById("city");
  const kEl = document.getElementById("country");

  function norm(x){ return (x||"").toString().toLowerCase(); }

  function render(){
    const q = norm(qEl.value);
    const st = norm(sEl.value);
    const city = norm(cEl.value);
    const country = norm(kEl.value);

    const filtered = records.filter(r => {
      const blob = norm((r.title||"") + " " + (r.paper||"") + " " + (r.sra_link||""));
      const okQ = !q || blob.includes(q);
      const okS = !st || norm(r.study_type) === st;
      const okC = !city || norm(r.city).includes(city);
      const okK = !country || norm(r.country).includes(country);
      return okQ && okS && okC && okK;
    });

    if (!filtered.length) {
      list.innerHTML = "<p class='muted'>No records match filters for this run.</p>";
      return;
    }

    list.innerHTML = filtered.slice().reverse().map(r => {
      const src = r.source || "";
      const st = r.study_type ? `<span class="tag">${r.study_type}</span>` : "";
      const loc = [r.city, r.country].filter(Boolean).join(", ");
      const locHtml = loc ? `<div class="muted" style="margin-top:6px">📍 ${loc}</div>` : "";
      const title = r.title ? `<div style="margin-top:8px">${r.title}</div>` : "";
      const paper = r.paper ? `<div><a href="${r.paper}" target="_blank"><b>Paper</b></a></div>` : "";
      const dataset = r.sra_link ? `<div><a href="${r.sra_link}" target="_blank"><b>Dataset</b></a></div>` : "";
      const uid = r.sra_uid ? `<div class="muted">SRA UID: <code>${r.sra_uid}</code></div>` : "";
      const pmid = r.pmid ? `<div class="muted">PMID: <code>${r.pmid}</code></div>` : "";
      return `<div class="card">
        <div class="muted">${src} ${st}</div>
        ${paper}
        ${dataset}
        ${uid}
        ${pmid}
        ${locHtml}
        ${title}
      </div>`;
    }).join("");
  }

  qEl.addEventListener("input", render);
  sEl.addEventListener("change", render);
  cEl.addEventListener("input", render);
  kEl.addEventListener("input", render);

  // If no records, still show a friendly message
  if (!records.length) {
    list.innerHTML = "<p class='muted'>No new linked datasets in this window.</p>";
  } else {
    render();
  }
})();
</script>
</body>
</html>
"""
    with open(DOCS_INDEX_HTML, "w", encoding="utf-8") as f:
        f.write(html)


def write_latest_json(records: List[Dict[str, Any]], query: str) -> None:
    payload = {
        "generated_at_utc": dt.datetime.utcnow().isoformat(timespec="seconds"),
        "query": query,
        "added_records": len(records),
        "records": records[-300:],  # keep Pages payload modest
    }
    with open(LATEST_JSON_PATH, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)


def write_latest_csv(records: List[Dict[str, Any]]) -> None:
    with open(LATEST_CSV_PATH, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        w.writeheader()
        for r in records[-2000:]:
            row = {k: r.get(k, "") for k in CSV_FIELDS}
            if isinstance(row.get("sra_uids"), list):
                row["sra_uids"] = " ".join(row["sra_uids"])
            w.writerow(row)


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

    p_daily = sub.add_parser("daily", help="Daily incremental update (last N days)")
    p_daily.add_argument("--days", type=int, default=1)
    p_daily.add_argument("--query", default=DEFAULT_QUERY)
    p_daily.add_argument("--max-per-day", type=int, default=200)

    args = ap.parse_args()
    ensure_dirs()

    seen = load_seen(SEEN_PATH)

    added_records_all: List[Dict[str, Any]] = []
    touched_years: Set[int] = set()
    added_ids_all: List[str] = []

    if args.cmd == "backfill-year":
        year = args.year
        start = dt.date(year, 1, 1)
        end = dt.date(year, 12, 31)
        cat_path = catalog_path_for_year(year)

        for day in daterange(start, end):
            recs, new_ids, tys = run_day(args.query, day, args.max_per_day, seen)
            if recs:
                append_jsonl(cat_path, recs)
                added_records_all.extend(recs)
            if new_ids:
                append_seen(SEEN_PATH, new_ids)
                for x in new_ids:
                    seen.add(x)
                added_ids_all.extend(new_ids)
            touched_years |= tys
            print(f"[{day.isoformat()}] +{len(recs)} records, +{len(new_ids)} ids")

    else:
        # daily
        end = dt.datetime.utcnow().date()
        start = end - dt.timedelta(days=max(args.days, 1))
        cat_path = catalog_path_for_year(end.year)

        for day in daterange(start, end):
            recs, new_ids, tys = run_day(args.query, day, args.max_per_day, seen)
            if recs:
                append_jsonl(cat_path, recs)
                added_records_all.extend(recs)
            if new_ids:
                append_seen(SEEN_PATH, new_ids)
                for x in new_ids:
                    seen.add(x)
                added_ids_all.extend(new_ids)
            touched_years |= tys
            print(f"[{day.isoformat()}] +{len(recs)} records, +{len(new_ids)} ids")

    # Write latest artifacts (Pages)
    write_latest_json(added_records_all, args.query)
    write_latest_csv(added_records_all)
    write_dashboard_html()

    # CSV exports (per-year) for any touched years (keeps it fast)
    for y in sorted(touched_years):
        out_csv = write_year_csv(y)
        if out_csv:
            print(f"Wrote {out_csv}")

    # Archive index + HTML
    idx = build_archive_index()
    write_archive_html(idx)

    print(
        f"UrbanScope done. added_records={len(added_records_all)} "
        f"added_ids={len(added_ids_all)} years_touched={sorted(touched_years) if touched_years else []}"
    )
    db_index = write_database_exports()
    print(f"DB exports updated: total_records={db_index['total_records']}")
  
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
