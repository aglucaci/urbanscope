#!/usr/bin/env python3
"""
urbanscope_radar.py

Create a clean, publication-ready radar chart for an "UrbanScope" style city/site profile.

What it does
- Loads per-city metrics from a CSV/TSV or JSON
- Selects one city (or the most recent / first row)
- Normalizes metrics (optional) and plots a radar/spider chart
- Saves PNG/SVG and (optionally) a JSON "card" with the values used

Input formats
1) CSV/TSV wide table (recommended):
   city,date,metric_a,metric_b,metric_c,...
   NYC,2026-01-21,0.31,120,0.82,...

2) JSON:
   {
     "city": "NYC",
     "date": "2026-01-21",
     "metrics": {"metric_a": 0.31, "metric_b": 120, "metric_c": 0.82}
   }

Normalization
- If you provide a ranges file, values are scaled to [0,1] using min/max.
- If no ranges file, you can still plot raw values (not ideal for radar).
- Ranges file is CSV/TSV with columns: metric,min,max
"""

from __future__ import annotations

import argparse
import json
import math
import os
from dataclasses import dataclass
from datetime import datetime
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


@dataclass
class RadarConfig:
    title: str
    label_wrap: int = 18
    fill_alpha: float = 0.18
    line_width: float = 2.5
    grid_alpha: float = 0.25
    tick_label_size: int = 9
    axis_label_size: int = 10
    title_size: int = 13
    dpi: int = 300


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate an UrbanScope radar chart from city metrics.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-i", "--input", required=True, help="Input CSV/TSV or JSON.")
    p.add_argument("--sep", default=None, help="CSV separator. Auto if omitted (comma/tsv).")
    p.add_argument("--city", default=None, help="City/site name to plot (matches 'city' column).")
    p.add_argument("--date", default=None, help="Optional date filter (matches 'date' column).")
    p.add_argument(
        "--metrics",
        default=None,
        help="Comma-separated list of metrics to include (in order). If omitted, uses all numeric columns.",
    )
    p.add_argument(
        "--rename",
        default=None,
        help="Optional JSON mapping to rename metrics for display, e.g. '{\"pm25\":\"PM2.5\"}'.",
    )
    p.add_argument(
        "--ranges",
        default=None,
        help="Optional ranges CSV/TSV with columns metric,min,max for normalization.",
    )
    p.add_argument(
        "--no-normalize",
        action="store_true",
        help="Disable normalization even if ranges are provided; plot raw values.",
    )
    p.add_argument(
        "--direction",
        default=None,
        help=(
            "Optional JSON mapping metric -> 'higher_better' or 'lower_better'. "
            "If lower_better, normalized value is inverted (1 - x)."
        ),
    )
    p.add_argument("-o", "--out", default="outputs", help="Output directory.")
    p.add_argument("--basename", default=None, help="Output base filename (without extension).")
    p.add_argument("--format", default="png", choices=["png", "svg", "pdf"], help="Figure format.")
    p.add_argument("--also-save-json", action="store_true", help="Also write a JSON card of values used.")
    p.add_argument("--title", default=None, help="Custom plot title.")
    p.add_argument("--show", action="store_true", help="Show plot interactively.")
    return p.parse_args()


def detect_sep(path: str, user_sep: Optional[str]) -> str:
    if user_sep is not None:
        return user_sep
    # simple heuristic: tsv if endswith .tsv or .tab
    lower = path.lower()
    if lower.endswith(".tsv") or lower.endswith(".tab"):
        return "\t"
    return ","


def wrap_label(s: str, width: int) -> str:
    if width <= 0 or len(s) <= width:
        return s
    parts = s.split(" ")
    if len(parts) == 1:
        # hard wrap
        return "\n".join([s[i : i + width] for i in range(0, len(s), width)])
    lines: List[str] = []
    cur: List[str] = []
    cur_len = 0
    for w in parts:
        if cur_len + len(w) + (1 if cur else 0) <= width:
            cur.append(w)
            cur_len += len(w) + (1 if cur_len else 0)
        else:
            lines.append(" ".join(cur))
            cur = [w]
            cur_len = len(w)
    if cur:
        lines.append(" ".join(cur))
    return "\n".join(lines)


def load_json_input(path: str) -> Tuple[str, Optional[str], Dict[str, float]]:
    with open(path, "r", encoding="utf-8") as f:
        obj = json.load(f)
    if "metrics" in obj and isinstance(obj["metrics"], dict):
        metrics = {k: float(v) for k, v in obj["metrics"].items()}
        city = str(obj.get("city", ""))
        date = obj.get("date", None)
        date = str(date) if date is not None else None
        return city, date, metrics

    # Allow a flat dict: {"city":..,"date":..,"metric_a":..}
    city = str(obj.get("city", ""))
    date = obj.get("date", None)
    date = str(date) if date is not None else None
    metrics = {}
    for k, v in obj.items():
        if k in {"city", "date"}:
            continue
        try:
            metrics[k] = float(v)
        except Exception:
            pass
    return city, date, metrics


def load_table_input(path: str, sep: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=sep)
    # normalize common column names
    cols = {c: c.strip() for c in df.columns}
    df.rename(columns=cols, inplace=True)
    return df


def pick_row(df: pd.DataFrame, city: Optional[str], date: Optional[str]) -> pd.Series:
    work = df.copy()
    if "city" in work.columns and city is not None:
        work = work[work["city"].astype(str) == str(city)]
    if "date" in work.columns and date is not None:
        work = work[work["date"].astype(str) == str(date)]
    if work.empty:
        raise ValueError("No rows matched your --city/--date filters.")
    # Prefer most recent date if available
    if "date" in work.columns:
        try:
            tmp = work.copy()
            tmp["_parsed_date"] = pd.to_datetime(tmp["date"], errors="coerce")
            tmp = tmp.sort_values("_parsed_date", ascending=False)
            return tmp.drop(columns=["_parsed_date"]).iloc[0]
        except Exception:
            pass
    return work.iloc[0]


def extract_metrics(
    row: pd.Series,
    metrics_list: Optional[List[str]],
) -> Dict[str, float]:
    if metrics_list:
        out = {}
        for m in metrics_list:
            if m not in row.index:
                raise KeyError(f"Metric '{m}' not found in input columns.")
            out[m] = float(row[m])
        return out

    # default: all numeric columns excluding city/date
    exclude = {"city", "date"}
    out = {}
    for k, v in row.items():
        if k in exclude:
            continue
        try:
            fv = float(v)
        except Exception:
            continue
        if math.isfinite(fv):
            out[str(k)] = fv
    if not out:
        raise ValueError("No numeric metrics found to plot.")
    return out


def load_ranges(path: str, sep: str) -> Dict[str, Tuple[float, float]]:
    rdf = pd.read_csv(path, sep=sep)
    needed = {"metric", "min", "max"}
    if not needed.issubset(set(rdf.columns)):
        raise ValueError(f"Ranges file must contain columns: {sorted(needed)}")
    ranges = {}
    for _, r in rdf.iterrows():
        m = str(r["metric"])
        mn = float(r["min"])
        mx = float(r["max"])
        if mx <= mn:
            raise ValueError(f"Invalid range for {m}: max <= min ({mn}, {mx})")
        ranges[m] = (mn, mx)
    return ranges


def normalize_metrics(
    metrics: Dict[str, float],
    ranges: Dict[str, Tuple[float, float]],
    direction: Optional[Dict[str, str]] = None,
) -> Dict[str, float]:
    norm = {}
    for m, v in metrics.items():
        if m not in ranges:
            raise KeyError(f"Metric '{m}' missing from ranges file.")
        mn, mx = ranges[m]
        x = (v - mn) / (mx - mn)
        x = float(np.clip(x, 0.0, 1.0))
        if direction and m in direction:
            if direction[m] == "lower_better":
                x = 1.0 - x
        norm[m] = x
    return norm


def radar_angles(n: int) -> np.ndarray:
    return np.linspace(0, 2 * np.pi, n, endpoint=False)


def make_radar_plot(
    labels: List[str],
    values: List[float],
    title: str,
    cfg: RadarConfig,
) -> plt.Figure:
    n = len(labels)
    ang = radar_angles(n)
    # close the loop
    ang_c = np.concatenate([ang, [ang[0]]])
    val_c = np.concatenate([np.array(values, dtype=float), [values[0]]])

    fig = plt.figure(figsize=(7.2, 7.2))
    ax = fig.add_subplot(111, polar=True)

    ax.set_theta_offset(np.pi / 2.0)
    ax.set_theta_direction(-1)

    ax.plot(ang_c, val_c, linewidth=cfg.line_width)
    ax.fill(ang_c, val_c, alpha=cfg.fill_alpha)

    ax.set_xticks(ang)
    ax.set_xticklabels(labels, fontsize=cfg.axis_label_size)

    # radial settings: assume normalized [0,1]
    ax.set_ylim(0, 1)
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(["0.2", "0.4", "0.6", "0.8", "1.0"], fontsize=cfg.tick_label_size)
    ax.yaxis.grid(True, alpha=cfg.grid_alpha)
    ax.xaxis.grid(True, alpha=cfg.grid_alpha)

    ax.set_title(title, fontsize=cfg.title_size, pad=18)
    fig.tight_layout()
    return fig


def safe_slug(s: str) -> str:
    s = s.strip().replace(" ", "_")
    return "".join(ch for ch in s if ch.isalnum() or ch in {"_", "-", "."})


def main() -> None:
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    in_lower = args.input.lower()
    rename_map = json.loads(args.rename) if args.rename else {}
    direction_map = json.loads(args.direction) if args.direction else None

    if args.metrics:
        metrics_list = [m.strip() for m in args.metrics.split(",") if m.strip()]
    else:
        metrics_list = None

    if in_lower.endswith(".json"):
        city0, date0, metrics0 = load_json_input(args.input)
        city = args.city or city0 or "Unknown"
        date = args.date or date0
        metrics = metrics0
    else:
        sep = detect_sep(args.input, args.sep)
        df = load_table_input(args.input, sep=sep)
        row = pick_row(df, args.city, args.date)
        city = str(row["city"]) if "city" in row.index else (args.city or "Unknown")
        date = str(row["date"]) if "date" in row.index else (args.date or None)
        metrics = extract_metrics(row, metrics_list)

    # order metrics if user specified
    if metrics_list:
        ordered_keys = metrics_list
    else:
        ordered_keys = list(metrics.keys())

    # normalization
    used_values = metrics.copy()
    normalized = False
    if args.ranges and not args.no_normalize:
        rsep = detect_sep(args.ranges, args.sep)
        ranges = load_ranges(args.ranges, sep=rsep)
        used_values = normalize_metrics(metrics, ranges, direction=direction_map)
        normalized = True

    # labels with optional rename + wrapping
    labels = []
    values = []
    for k in ordered_keys:
        disp = rename_map.get(k, k)
        disp = wrap_label(str(disp), width=RadarConfig(title="").label_wrap)
        labels.append(disp)
        values.append(float(used_values[k]))

    # Title + basename
    stamp = datetime.now().strftime("%Y-%m-%d")
    date_part = date if date else stamp
    title = args.title or f"UrbanScope Radar — {city} ({date_part})"
    base = args.basename or f"urbanscope_radar_{safe_slug(city)}_{safe_slug(date_part)}"

    cfg = RadarConfig(title=title)
    fig = make_radar_plot(labels, values, title=title, cfg=cfg)

    out_path = os.path.join(args.out, f"{base}.{args.format}")
    fig.savefig(out_path, dpi=cfg.dpi, bbox_inches="tight")
    plt.close(fig)

    if args.also_save_json:
        card = {
            "city": city,
            "date": date_part,
            "normalized": normalized,
            "metrics": {k: float(used_values[k]) for k in ordered_keys},
            "raw_metrics": {k: float(metrics[k]) for k in ordered_keys},
        }
        card_path = os.path.join(args.out, f"{base}.json")
        with open(card_path, "w", encoding="utf-8") as f:
            json.dump(card, f, indent=2)

    if args.show:
        # reload and display quickly
        img = plt.imread(out_path) if args.format == "png" else None
        if img is not None:
            plt.figure()
            plt.imshow(img)
            plt.axis("off")
            plt.show()

    print(out_path)


if __name__ == "__main__":
    main()
