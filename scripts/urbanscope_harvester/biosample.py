from __future__ import annotations
import re
import xml.etree.ElementTree as ET
from typing import Any, Dict, List

from .ncbi import http_get, _eutils_params
from .config import EUTILS
from .utils import _norm

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

def efetch_biosample_xml(accession_or_uid: str):
    params = _eutils_params({"db": "biosample", "id": accession_or_uid, "retmode": "xml"})
    url = EUTILS + "efetch.fcgi?" + __import__("urllib.parse").parse.urlencode(params)
    xmltxt = http_get(url).decode(errors="replace")
    return xmltxt, url

def parse_biosample_attributes_from_xml(xmltxt: str) -> Dict[str, Any]:
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

def biosample_display_card(bsd: Dict[str, Any]) -> Dict[str, Any]:
    """
    Make a small, UI-friendly subset for HTML display.
    Keeps the raw attributes, but also surfaces a few common ones.
    """
    if not isinstance(bsd, dict) or not bsd:
        return {}

    attrs = bsd.get("attributes", {}) if isinstance(bsd.get("attributes", {}), dict) else {}

    def pick(*keys):
        for k in keys:
            if k in attrs and str(attrs[k]).strip():
                return str(attrs[k]).strip()
        return ""

    # Common keys seen across BioSample records
    env_biome = pick("env_biome", "environment (biome)", "environment_biome")
    env_feature = pick("env_feature", "environment (feature)", "environment_feature")
    env_material = pick("env_material", "environment (material)", "environment_material")
    host = pick("host", "host scientific name", "host_taxid", "host_common_name")
    samp_type = pick("sample_type", "sample type", "isolation_source", "isolation source")
    collection_date = pick("collection_date", "collection date")
    sample_name = pick("sample_name", "sample name", "sample_name_alias")
    depth = pick("depth", "elevation", "altitude")
    temp = pick("temperature", "temp")
    ph = pick("ph")

    return {
        "accession": bsd.get("accession", ""),
        "title": bsd.get("title", ""),
        "organism": bsd.get("organism", ""),
        "collection_date": collection_date,
        "sample_name": sample_name,
        "sample_type": samp_type,
        "host": host,
        "env_biome": env_biome,
        "env_feature": env_feature,
        "env_material": env_material,
        "depth_or_altitude": depth,
        "temperature": temp,
        "ph": ph,
        "attributes": attrs,         # keep full attrs for deep UI drill-down
        "efetch_url": bsd.get("efetch_url", ""),
        "error": bsd.get("error", ""),
    }


def infer_geo(biosample_details: Dict[str, Any], fallbacks: List[str]) -> Dict[str, Any]:
    attrs = (biosample_details or {}).get("attributes", {}) if isinstance(biosample_details, dict) else {}
    raw_geo = ""

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
