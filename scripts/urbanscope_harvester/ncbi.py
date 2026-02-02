from __future__ import annotations
import urllib.request, urllib.parse
import xml.etree.ElementTree as ET
from typing import Any, Dict, List, Tuple

from .config import EUTILS, NCBI_API_KEY, TOOL_NAME, NCBI_EMAIL, BIOPROJECT_RE
from .utils import _sleep_backoff

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

def _eutils_params(extra: Dict[str, str]) -> Dict[str, str]:
    p = dict(extra)
    if NCBI_API_KEY:
        p["api_key"] = NCBI_API_KEY
    if TOOL_NAME:
        p["tool"] = TOOL_NAME
    if NCBI_EMAIL:
        p["email"] = NCBI_EMAIL
    return p

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

def esearch_day(db: str, term: str, day: str, retmax: int, datetype: str = "edat") -> Tuple[List[str], str]:
    params = _eutils_params({
        "db": db, "term": term, "retmode": "xml",
        "mindate": day, "maxdate": day, "datetype": datetype,
        "retmax": str(retmax),
    })
    url = EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(params)
    root = parse_xml(http_get(url))
    ids = [x.text for x in root.findall(".//Id") if x.text]
    return ids, url

def esearch_recent(db: str, term: str, reldate_days: int, retmax: int, datetype: str = "edat") -> Tuple[List[str], str]:
    params = _eutils_params({
        "db": db, "term": term, "retmode": "xml",
        "reldate": str(reldate_days), "datetype": datetype,
        "retmax": str(retmax), "sort": "date",
    })
    url = EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(params)
    root = parse_xml(http_get(url))
    ids = [x.text for x in root.findall(".//Id") if x.text]
    return ids, url

def esearch_history(db: str, term: str, retstart: int, retmax: int, sort: str = ""):
    params = _eutils_params({
        "db": db, "term": term, "retmode": "xml",
        "retstart": str(retstart), "retmax": str(retmax),
        "usehistory": "n",
    })
    if sort:
        params["sort"] = sort
    url = EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(params)
    root = parse_xml(http_get(url))
    ids = [x.text for x in root.findall(".//Id") if x.text]
    count_total = int((root.findtext(".//Count") or "0").strip() or "0")
    return ids, count_total, url

def efetch_runinfo_text(uid: str) -> Tuple[str, str]:
    params = _eutils_params({"db": "sra", "id": uid, "rettype": "runinfo", "retmode": "text"})
    url = EUTILS + "efetch.fcgi?" + urllib.parse.urlencode(params)
    text = http_get(url).decode(errors="replace")
    return text, url

def esummary_sra(uids: List[str]) -> Tuple[Dict[str, Dict[str, Any]], str]:
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
