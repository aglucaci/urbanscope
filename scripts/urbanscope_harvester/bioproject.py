from __future__ import annotations
from typing import Any, Dict, Tuple, Optional

from .ncbi import esearch_any, esummary
from .config import BIOPROJECT_RE

from xml.etree import ElementTree as ET
from pathlib import Path

def _txt(node, tag: str, default: str = "") -> str:
    x = node.find(tag) if node is not None else None
    return (x.text or "").strip() if x is not None and x.text is not None else default

def _attr(node, path: str, attr: str, default: str = "") -> str:
    x = node.find(path) if node is not None else None
    return (x.get(attr) or "").strip() if x is not None else default

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
    root, _ = esummary("bioproject", [uid])

    # -------------------------------
    # FORMAT 1: Rich DocumentSummary
    # -------------------------------
    doc = root.find(".//DocumentSummary/Project")
    if doc is not None:
        ds = root.find(".//DocumentSummary")

        acc = ds.find("./Project/ProjectID/ArchiveID").get("accession", "").strip().upper()
        title = _txt(ds, "./Project/ProjectDescr/Title")
        desc = _txt(ds, "./Project/ProjectDescr/Description")

        data_type = _txt(
            ds,
            "./Project/ProjectType/ProjectTypeSubmission/IntendedDataTypeSet/DataType",
        )
        if not data_type:
            data_type = ds.find(
                "./Project/ProjectType/ProjectTypeSubmission/Objectives/Data"
            ).get("data_type", "").strip()

        submission_date = ds.find("./Submission").get("submitted", "").strip()
        last_update = ds.find("./Submission").get("last_update", "").strip()
        center = _txt(ds, "./Submission/Description/Organization/Name")

        return {
            "uid": uid,
            "accession": acc,
            "title": title,
            "description": desc,
            "organism": "",
            "data_type": data_type,
            "submission_date": submission_date,
            "last_update": last_update,
            "center_name": center,
            "ncbi": {
                "bioproject_uid": uid,
                "bioproject_url": f"https://www.ncbi.nlm.nih.gov/bioproject/{uid}",
            },
        }

    # ------------------------------------------------
    # FORMAT 2: Flat DocumentSummary (your example)
    # ------------------------------------------------
    doc = root.find(".//DocumentSummary")
    if doc is not None and doc.find("Project_Acc") is not None:
        acc = _txt(doc, "Project_Acc").upper()
        title = _txt(doc, "Project_Title")
        desc = _txt(doc, "Project_Description")
        organism = _txt(doc, "Organism_Name")
        data_type = _txt(doc, "Project_Data_Type")
        submission_date = _txt(doc, "Registration_Date")
        last_update = ""  # not exposed in this variant

        # Prefer primary submitter org
        center = _txt(doc, "Submitter_Organization")
        if not center:
            orgs = doc.find("Submitter_Organization_List")
            if orgs is not None:
                vals = [x.text.strip() for x in orgs.findall("string") if x.text]
                center = vals[0] if vals else ""

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
        }

    # --------------------------------
    # FORMAT 3: Legacy DocSum / Item
    # --------------------------------
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



def get_bioproject_details(accession: str, bp_cache: Dict[str, Any], uid_cache: Dict[str, str]) -> Dict[str, Any]:
    accession = (accession or "").strip().upper()
    if not accession:
        return {}
    if accession in bp_cache:
        return bp_cache[accession] or {}

    if not BIOPROJECT_RE.match(accession):
        bp_cache[accession] = {"accession": accession, "uid": "", "error": "invalid_accession"}
        return bp_cache[accession]

    uid = bioproject_accession_to_uid(accession, uid_cache)
    if not uid:
        bp_cache[accession] = {"accession": accession, "uid": "", "error": "uid_not_found"}
        return bp_cache[accession]

    print("[INFO] Parsing Bioproject esummary:", uid)
    details = parse_bioproject_esummary(uid)
    details["accession"] = details.get("accession") or accession
    bp_cache[accession] = details
    print("details:", details)
    return details
