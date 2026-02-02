from __future__ import annotations
import os, re, glob
from typing import Any, Dict, Iterator, List

from .config import DATA_DIR, DB_DIR, DOCS_LATEST_SRR, MAX_OUTPUT_BYTES
from .utils import (
    utc_now, read_json, write_json, iter_jsonl_glob,
    write_json_array_chunked, file_size
)

def _find_year_catalog_prefixes() -> List[str]:
    prefixes = set()
    for p in glob.glob(os.path.join(DATA_DIR, "srr_catalog_*.jsonl")):
        prefixes.add(p)
    for p in glob.glob(os.path.join(DATA_DIR, "srr_catalog_*_part[0-9][0-9][0-9].jsonl")):
        base = re.sub(r"_part[0-9]{3}\.jsonl$", ".jsonl", p)
        prefixes.add(base)

    def year_key(path: str) -> int:
        m = re.search(r"srr_catalog_(\d{4})", os.path.basename(path))
        return int(m.group(1)) if m else 0

    return sorted(prefixes, key=year_key)

def rebuild_srr_exports_chunked():
    prefixes = _find_year_catalog_prefixes()

    def all_records() -> Iterator[Dict[str, Any]]:
        for base in prefixes:
            yield from iter_jsonl_glob(base)

    manifest = write_json_array_chunked(
        out_prefix=os.path.join(DB_DIR, "srr_records"),
        records_iter=all_records(),
        max_bytes=MAX_OUTPUT_BYTES,
    )

    years = set()
    for p in prefixes:
        m = re.search(r"srr_catalog_(\d{4})", os.path.basename(p))
        if m:
            years.add(int(m.group(1)))
    manifest["years"] = sorted(years)

    write_json(os.path.join(DB_DIR, "srr_records_manifest.json"), manifest)
    write_json(os.path.join(DB_DIR, "srr_index.json"), {
        "generated": utc_now(),
        "total_srr_records": manifest["total_records"],
        "parts": [os.path.basename(x["path"]) for x in manifest["parts"]],
        "years": manifest["years"],
    })

    from .config import BIOPROJECT_CACHE, BIOSAMPLE_CACHE
    bp_cache = read_json(BIOPROJECT_CACHE, {})
    if bp_cache:
        write_json(os.path.join(DB_DIR, "bioprojects.json"), bp_cache)
    bs_cache = read_json(BIOSAMPLE_CACHE, {})
    if bs_cache:
        write_json(os.path.join(DB_DIR, "biosamples.json"), bs_cache)

def write_latest_srr_safe(latest_items: List[Dict[str, Any]]):
    payload = {"generated_utc": utc_now(), "count": len(latest_items), "items": latest_items}

    tmp = DOCS_LATEST_SRR + ".tmp"
    import json
    os.makedirs(os.path.dirname(DOCS_LATEST_SRR) or ".", exist_ok=True)
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)
    sz = file_size(tmp)

    if sz <= MAX_OUTPUT_BYTES:
        os.replace(tmp, DOCS_LATEST_SRR)
        return

    lo, hi = 0, len(latest_items)
    best = 0
    while lo <= hi:
        mid = (lo + hi) // 2
        test = {"generated_utc": utc_now(), "count": len(latest_items), "items": latest_items[:mid]}
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(test, f, ensure_ascii=False, indent=2)
        if file_size(tmp) <= MAX_OUTPUT_BYTES:
            best = mid
            lo = mid + 1
        else:
            hi = mid - 1

    final = {"generated_utc": utc_now(), "count": len(latest_items), "items": latest_items[:best]}
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(final, f, ensure_ascii=False, indent=2)
    os.replace(tmp, DOCS_LATEST_SRR)
