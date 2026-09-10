"""Compare two sealed runs without printing or publishing scientific results."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest().upper()


def numeric_leaves(value, prefix=""):
    out = {}
    if isinstance(value, dict):
        for key, item in value.items():
            out.update(numeric_leaves(item, f"{prefix}.{key}" if prefix else key))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            out.update(numeric_leaves(item, f"{prefix}[{index}]"))
    elif isinstance(value, (int, float)) and not isinstance(value, bool):
        out[prefix] = float(value)
    return out


def csv_numeric(path: Path):
    values = {}
    with path.open(newline="", encoding="utf-8-sig") as f:
        for row_index, row in enumerate(csv.DictReader(f)):
            for key, raw in row.items():
                try:
                    values[f"{row_index}.{key}"] = float(raw)
                except (TypeError, ValueError):
                    pass
    return values


def max_diff(a, b):
    if set(a) != set(b):
        return None
    return max((abs(a[key] - b[key]) for key in a), default=0.0)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run-a", type=Path, required=True)
    p.add_argument("--run-b", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    files = [
        "data/restricted_nigeria.csv", "data/restricted_burkina_faso.csv", "data/derivation_manifest.json",
        "sealed/sealed_results.json", "sealed/sealed_coefficients.csv", "sealed/sealed_covariances.csv",
        "sealed/sealed_estimates.csv", "sealed/R_sessionInfo.txt",
    ]
    hashes = {}
    for name in files:
        ha, hb = sha256(args.run_a / name), sha256(args.run_b / name)
        hashes[name] = {"sha256": ha, "exact_match": ha == hb}
    ja = json.loads((args.run_a / "sealed/sealed_results.json").read_text(encoding="utf-8"))
    jb = json.loads((args.run_b / "sealed/sealed_results.json").read_text(encoding="utf-8"))
    numeric = numeric_leaves(ja), numeric_leaves(jb)
    coefficient_diff = max_diff(csv_numeric(args.run_a / "sealed/sealed_coefficients.csv"), csv_numeric(args.run_b / "sealed/sealed_coefficients.csv"))
    covariance_diff = max_diff(csv_numeric(args.run_a / "sealed/sealed_covariances.csv"), csv_numeric(args.run_b / "sealed/sealed_covariances.csv"))
    estimate_diff = max_diff(csv_numeric(args.run_a / "sealed/sealed_estimates.csv"), csv_numeric(args.run_b / "sealed/sealed_estimates.csv"))
    payload = {
        "comparison_scope": ["analysis samples", "survey design metadata", "model coefficients", "coefficient covariance matrices", "standardized risks", "contrasts", "confidence intervals", "P values", "tables", "session information"],
        "criterion": "Exact SHA-256 equality is required; numeric maximum absolute differences are also computed and must be zero.",
        "artifact_hashes": hashes,
        "numeric_max_absolute_difference": max_diff(*numeric),
        "coefficient_max_absolute_difference": coefficient_diff,
        "covariance_max_absolute_difference": covariance_diff,
        "estimate_table_max_absolute_difference": estimate_diff,
        "pass": all(item["exact_match"] for item in hashes.values()) and max_diff(*numeric) == 0 and coefficient_diff == 0 and covariance_diff == 0 and estimate_diff == 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
