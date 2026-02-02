from __future__ import annotations
from typing import Any, Dict, List
from .utils import _norm

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

    if "amplicon" in strat or "amplicon" in blob:
        hits.append("amplicon"); tags.append("amplicon")
        if "16s" in blob or ("rrna" in blob and "16s" in blob):
            hits.append("16s"); tags.append("16S")
            return {"assay_class": "16S", "assay_tags": tags, "confidence": "high", "rationale": hits}
        if "its" in blob:
            hits.append("its"); tags.append("ITS")
            return {"assay_class": "ITS", "assay_tags": tags, "confidence": "high", "rationale": hits}
        return {"assay_class": "Amplicon", "assay_tags": tags, "confidence": "high", "rationale": hits}

    if strat in ("rna-seq", "transcriptome") or "rna-seq" in blob or "metatranscriptom" in blob:
        hits.append("rna-seq/metatranscriptome"); tags.append("RNA")
        return {"assay_class": "RNA-seq", "assay_tags": tags, "confidence": "high", "rationale": hits}

    if strat in ("wgs", "metagenomic") or "shotgun" in blob or "wgs" in blob or "metagenom" in blob:
        hits.append("wgs/shotgun/metagenomic"); tags.append("shotgun")
        return {"assay_class": "WGS", "assay_tags": tags, "confidence": "high", "rationale": hits}

    if "pcr" in sel or "rrna" in sel:
        hits.append("PCR/rRNA selection"); tags.append("targeted")
        if "16s" in blob:
            hits.append("16s"); tags.append("16S")
            return {"assay_class": "16S", "assay_tags": tags, "confidence": "medium", "rationale": hits}
        if "its" in blob:
            hits.append("its"); tags.append("ITS")
            return {"assay_class": "ITS", "assay_tags": tags, "confidence": "medium", "rationale": hits}
        return {"assay_class": "Amplicon", "assay_tags": tags, "confidence": "medium", "rationale": hits}

    return {"assay_class": "Unknown", "assay_tags": [], "confidence": "low", "rationale": []}
