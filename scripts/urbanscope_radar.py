#!/usr/bin/env python3
"""
UrbanScope — Urban Metagenomics & Metatranscriptomics Dataset Radar

- Finds recent urban/built-environment metagenomics studies
- Links PubMed papers to SRA datasets when available
- 100% free: NCBI E-utilities only
- Output designed for GitHub Releases or Pages

Usage:
  python urban_dataset_radar.py --days 1 --max 50
"""

import datetime as dt
import json
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from typing import List, Dict

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

QUERY = (
    '(urban OR city OR cities OR subway OR transit OR "built environment" '
    'OR wastewater OR sewage) AND '
    '(metagenomics OR metatranscriptomics)'
)

def get(url: str) -> bytes:
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "urbanscope/1.0"}
    )
    with urllib.request.urlopen(req, timeout=30) as r:
        return r.read()

def esearch(db: str, term: str, mindate: str, maxdate: str, retmax: int) -> List[str]:
    params = {
        "db": db,
        "term": term,
        "retmode": "xml",
        "retmax": str(retmax),
        "mindate": mindate,
        "maxdate": maxdate,
        "datetype": "pdat",
    }
    url = EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(params)
    root = ET.fromstring(get(url))
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
    root = ET.fromstring(get(url))

    links = {}
    for linkset in root.findall("LinkSet"):
        pmid = linkset.findtext("IdList/Id")
        sra_ids = [
            l.text for l in linkset.findall(".//Link/Id") if l.text
        ]
        if pmid and sra_ids:
            links[pmid] = sra_ids
    return links

def main():
    today = dt.date.today()
    start = today - dt.timedelta(days=1)

    pmids = esearch(
        "pubmed",
        QUERY,
        start.strftime("%Y/%m/%d"),
        today.strftime("%Y/%m/%d"),
        retmax=50,
    )
    time.sleep(0.4)

    pmid_to_sra = elink_pubmed_to_sra(pmids)

    results = []
    for pmid, sra_ids in pmid_to_sra.items():
        results.append({
            "pmid": pmid,
            "paper": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            "sra_ids": sra_ids,
            "sra_link": f"https://www.ncbi.nlm.nih.gov/sra/?term={' OR '.join(sra_ids)}"
        })

    payload = {
        "generated_at": dt.datetime.utcnow().isoformat(),
        "query": QUERY,
        "count": len(results),
        "results": results,
    }

    with open("latest.json", "w") as f:
        json.dump(payload, f, indent=2)

    print(f"UrbanScope: found {len(results)} urban omics datasets")

if __name__ == "__main__":
    main()
