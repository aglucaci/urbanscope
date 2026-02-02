from __future__ import annotations
import os, json, time, random, datetime as dt, re
from typing import Any, Dict, Iterable, Iterator, List, Set

from .config import (
    DATA_DIR, DOCS_DIR, DB_DIR, CACHE_DIR, DEBUG_DIR, DOCS_DEBUG_DIR, MAX_OUTPUT_BYTES
)

def ensure_dirs():
    for p in [DATA_DIR, DOCS_DIR, DB_DIR, CACHE_DIR, DEBUG_DIR, DOCS_DEBUG_DIR]:
        os.makedirs(p, exist_ok=True)

def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")

def _sleep_backoff(i: int):
    time.sleep(0.6 * (2 ** i) + random.random() * 0.25)

def parse_int(x: str, default: int = 0) -> int:
    try:
        return int(x)
    except Exception:
        return default

def load_set(path: str) -> Set[str]:
    if not os.path.exists(path):
        return set()
    with open(path, encoding="utf-8") as f:
        return set(x.strip() for x in f if x.strip())

def append_lines(path: str, vals: Iterable[str]):
    vals = [v for v in vals if v]
    if not vals:
        return
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "a", encoding="utf-8") as f:
        for v in vals:
            f.write(v + "\n")

def read_json(path: str, default):
    if not os.path.exists(path):
        return default
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return default

def write_json(path: str, obj: Any):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)
    os.replace(tmp, path)

def _norm(s: str) -> str:
    return re.sub(r"\s+", " ", (s or "").strip().lower())

def inc(d: Dict[str, int], k: str, n: int = 1):
    d[k] = d.get(k, 0) + n

# ----------------------------
# Size-safe output helpers
# ----------------------------

def file_size(path: str) -> int:
    try:
        return os.path.getsize(path)
    except Exception:
        return 0

def rotating_path(base_path: str, max_bytes: int = MAX_OUTPUT_BYTES) -> str:
    """
    For JSONL appends. If base_path exceeds max_bytes, write to base_path_partNNN.jsonl.
    If base_path already has _partNNN, keep writing to that until full.
    """
    root, ext = os.path.splitext(base_path)
    if ext.lower() != ".jsonl":
        return base_path

    if not os.path.exists(base_path):
        return base_path

    if file_size(base_path) < max_bytes:
        return base_path

    i = 0
    while True:
        p = f"{root}_part{i:03d}{ext}"
        if (not os.path.exists(p)) or file_size(p) < max_bytes:
            return p
        i += 1

def append_jsonl(path: str, records: List[Dict[str, Any]], max_bytes: int = MAX_OUTPUT_BYTES):
    """
    Append records to JSONL, rotating files to keep each <= max_bytes.
    Rotation boundary is checked before each record write (safe even for big records).
    """
    if not records:
        return

    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    cur_path = rotating_path(path, max_bytes=max_bytes)
    f = open(cur_path, "a", encoding="utf-8")
    try:
        for r in records:
            line = json.dumps(r, ensure_ascii=False) + "\n"
            if file_size(cur_path) + len(line.encode("utf-8")) > max_bytes:
                f.close()
                cur_path = rotating_path(path, max_bytes=max_bytes)
                f = open(cur_path, "a", encoding="utf-8")
            f.write(line)
    finally:
        f.close()

def append_jsonl_one(path: str, rec: Dict[str, Any], max_bytes: int = MAX_OUTPUT_BYTES):
    append_jsonl(path, [rec], max_bytes=max_bytes)

def iter_jsonl(path: str) -> Iterator[Dict[str, Any]]:
    if not os.path.exists(path):
        return
    with open(path, encoding="utf-8") as fh:
        for ln in fh:
            ln = ln.strip()
            if ln:
                yield json.loads(ln)

def iter_jsonl_glob(prefix_path: str) -> Iterator[Dict[str, Any]]:
    """
    Iterate base.jsonl plus any base_partNNN.jsonl in order.
    """
    import glob
    root, ext = os.path.splitext(prefix_path)
    paths = []
    if os.path.exists(prefix_path):
        paths.append(prefix_path)
    paths.extend(sorted(glob.glob(f"{root}_part[0-9][0-9][0-9]{ext}")))
    for p in paths:
        yield from iter_jsonl(p)

def write_json_array_chunked(
    out_prefix: str,
    records_iter: Iterable[Dict[str, Any]],
    max_bytes: int = MAX_OUTPUT_BYTES,
) -> Dict[str, Any]:
    """
    Write JSON arrays into part files to keep each <= max_bytes.
    Produces out_prefix_part000.json, out_prefix_part001.json, ...
    Returns a manifest dict.
    """
    parts = []
    part_idx = 0

    def part_path(i: int) -> str:
        return f"{out_prefix}_part{i:03d}.json"

    os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)

    cur_path = part_path(part_idx)
    cur = open(cur_path, "w", encoding="utf-8")
    cur.write("[\n")
    first = True
    n_total = 0
    n_part = 0

    def cur_bytes() -> int:
        try:
            return cur.tell()
        except Exception:
            return file_size(cur_path)

    try:
        for rec in records_iter:
            blob = json.dumps(rec, ensure_ascii=False, indent=2)
            entry = ("" if first else ",\n") + blob
            if cur_bytes() + len(entry.encode("utf-8")) + len("\n]\n".encode("utf-8")) > max_bytes and not first:
                cur.write("\n]\n")
                cur.close()
                parts.append({"path": cur_path, "records": n_part})
                part_idx += 1
                cur_path = part_path(part_idx)
                cur = open(cur_path, "w", encoding="utf-8")
                cur.write("[\n")
                first = True
                n_part = 0
                entry = blob  # first entry

            cur.write(entry)
            first = False
            n_total += 1
            n_part += 1

        cur.write("\n]\n")
        cur.close()
        parts.append({"path": cur_path, "records": n_part})
    finally:
        try:
            cur.close()
        except Exception:
            pass

    return {"generated_utc": utc_now(), "total_records": n_total, "parts": parts}
