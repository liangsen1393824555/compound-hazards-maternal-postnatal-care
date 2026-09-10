"""Frozen H/D/K construction and aggregate support gate.

This is deliberately a descriptive proof-of-data step. It reports only
country-level counts by predeclared H/D/K groups; it creates no effect estimate,
association, model, test, P value, or public spatial microdata.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import tempfile
from datetime import date, datetime, timezone
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
import shapefile
from scipy.stats import gamma, norm

from internal_variable_denominator_gate import BASE, latest_most_recent_birth, mother_pnc_within_2_days


COUNTRIES = {"nigeria": ("2021-12", "2024-04"), "burkina_faso": ("2019-08", "2021-11")}
# The 50-km definition is primary.  The 25/100-km variants are predeclared
# sensitivities but deliberately not expanded during this proof-of-data phase.
RADII = (50,)


def atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
        temporary = Path(handle.name)
    os.replace(temporary, path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def add_month(value: str, delta: int) -> str:
    year, month = map(int, value.split("-"))
    ordinal = year * 12 + (month - 1) + delta
    return f"{ordinal // 12:04d}-{ordinal % 12 + 1:02d}"


def cmc_to_month(cmc: float) -> str:
    if not np.isfinite(cmc):
        raise ValueError("Missing outcome CMC")
    ordinal = int(cmc) - 1
    return f"{1900 + ordinal // 12:04d}-{ordinal % 12 + 1:02d}"


def month_start(value: str) -> date:
    year, month = map(int, value.split("-"))
    return date(year, month, 1)


def month_end(value: str) -> date:
    return month_start(add_month(value, 1)) - pd.Timedelta(days=1)


def cluster_key(value: object) -> str:
    return str(int(float(value)))


def protected_ge(work: Path, country: str) -> tuple[np.ndarray, dict[str, int]]:
    shp = next((work / "staging" / "official_dhs_members" / country / "ge").glob("*.shp"))
    reader = shapefile.Reader(str(shp))
    points: list[tuple[float, float]] = []
    mapping: dict[str, int] = {}
    for index, item in enumerate(reader.iterShapeRecords()):
        record = item.record.as_dict()
        key = cluster_key(record["DHSCLUST"])
        if key in mapping:
            raise RuntimeError(f"Duplicate protected GE cluster join key for {country}")
        mapping[key] = index
        points.append(tuple(item.shape.points[0]))
    if not points:
        raise RuntimeError(f"No protected GE points for {country}")
    return np.asarray(points, dtype="float64"), mapping


def load_npz(path: Path, value: str) -> tuple[list[str], np.ndarray]:
    with np.load(path) as source:
        months = [str(item) for item in source["months"].tolist()]
        matrix = np.asarray(source[value], dtype="float64")
    if matrix.ndim != 2 or matrix.shape[1] != len(months) or not np.isfinite(matrix).all():
        raise RuntimeError(f"Invalid protected value artifact: {path}")
    return months, matrix


def spi3(precip: np.ndarray, months: list[str]) -> tuple[np.ndarray, dict]:
    """SPI-3: rolling three-month total; month-specific 1991--2020 gamma CDF."""
    if len(months) != precip.shape[1]:
        raise RuntimeError("CHIRPS period/value matrix mismatch")
    totals = np.full_like(precip, np.nan, dtype="float64")
    totals[:, 2:] = precip[:, :-2] + precip[:, 1:-1] + precip[:, 2:]
    output = np.full_like(totals, np.nan, dtype="float64")
    baseline_n: list[int] = []
    for month_number in range(1, 13):
        all_positions = np.asarray([i for i, value in enumerate(months) if int(value[-2:]) == month_number], dtype=int)
        reference = np.asarray([i for i in all_positions if "1991-01" <= months[i] <= "2020-12" and np.isfinite(totals[0, i])], dtype=int)
        if len(reference) < 29:
            raise RuntimeError(f"Insufficient SPI baseline months for calendar month {month_number}")
        baseline_n.append(int(len(reference)))
        for row in range(precip.shape[0]):
            base = totals[row, reference]
            if not np.isfinite(base).all():
                raise RuntimeError("Invalid SPI baseline total")
            zero_fraction = float(np.mean(base == 0))
            positive = base[base > 0]
            if len(positive) < 3:
                raise RuntimeError("Insufficient positive baseline values for gamma SPI")
            shape, _, scale = gamma.fit(positive, floc=0)
            target = totals[row, all_positions]
            cdf = np.full(len(target), np.nan, dtype="float64")
            valid = np.isfinite(target)
            cdf[valid & (target <= 0)] = zero_fraction
            positive_target = valid & (target > 0)
            cdf[positive_target] = zero_fraction + (1 - zero_fraction) * gamma.cdf(target[positive_target], shape, loc=0, scale=scale)
            cdf[valid] = np.clip(cdf[valid], 1e-7, 1 - 1e-7)
            output[row, all_positions] = norm.ppf(cdf)
    return output, {
        "indicator": "SPI-3",
        "accumulation": "rolling three-month precipitation total ending in each calendar month",
        "baseline": "1991-01 through 2020-12; month-specific samples",
        "distribution": "gamma MLE on positive baseline totals with empirical zero probability q; CDF=q+(1-q)*GammaCDF for positive values",
        "baseline_sample_months_by_calendar_month": baseline_n,
        "threshold": "SPI-3 <= -1.0",
        "note": "The first two 1991 rolling totals are unavailable because retrieval begins in 1991-01; January and February therefore have 29 valid baseline endpoints, other calendar months have 30.",
    }


def heat(t2m: np.ndarray, months: list[str]) -> tuple[np.ndarray, dict]:
    reference = np.asarray([i for i, value in enumerate(months) if "1991-01" <= value <= "2020-12"], dtype=int)
    if len(reference) != 360:
        raise RuntimeError("ERA5-Land 1991--2020 baseline is incomplete")
    thresholds = np.empty((t2m.shape[0], 12), dtype="float64")
    result = np.zeros(t2m.shape, dtype=bool)
    for month_number in range(1, 13):
        cols = np.asarray([i for i in reference if int(months[i][-2:]) == month_number], dtype=int)
        thresholds[:, month_number - 1] = np.quantile(t2m[:, cols], 0.90, axis=1, method="linear")
        every = np.asarray([i for i, value in enumerate(months) if int(value[-2:]) == month_number], dtype=int)
        result[:, every] = t2m[:, every] > thresholds[:, month_number - 1, None]
    return result, {"indicator": "H", "definition": "At least one T2m month in the prior three full calendar months is strictly greater than its protected-gridpoint/calendar-month 1991--2020 P90.", "source_units": "K", "baseline": "1991-01 through 2020-12", "percentile": 90, "comparison": ">"}


def distance_km(points: np.ndarray, events: pd.DataFrame) -> np.ndarray:
    lat1 = np.deg2rad(points[:, 1])[:, None]
    lon1 = np.deg2rad(points[:, 0])[:, None]
    lat2 = np.deg2rad(events["latitude"].to_numpy(dtype="float64"))[None, :]
    lon2 = np.deg2rad(events["longitude"].to_numpy(dtype="float64"))[None, :]
    hav = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 6371.0088 * 2 * np.arcsin(np.minimum(1, np.sqrt(hav)))


def k_by_radius(points: np.ndarray, events: pd.DataFrame, outcome_months: list[str]) -> tuple[dict[int, np.ndarray], dict[int, int]]:
    distances = distance_km(points, events)
    value = {radius: np.zeros((len(points), len(outcome_months)), dtype=bool) for radius in RADII}
    supporting_events = {radius: np.zeros(len(events), dtype=bool) for radius in RADII}
    for col, outcome in enumerate(outcome_months):
        start, end = month_start(add_month(outcome, -3)), month_end(add_month(outcome, -1))
        active = (events["date_start"].dt.date <= end) & (events["date_end"].dt.date >= start)
        if not active.any():
            continue
        active_index = np.flatnonzero(active.to_numpy())
        for radius in RADII:
            near = distances[:, active_index] <= radius
            value[radius][:, col] = near.any(axis=1)
            supporting_events[radius][active_index] |= near.any(axis=0)
    return value, {radius: int(mask.sum()) for radius, mask in supporting_events.items()}


def group_table(frame: pd.DataFrame, radius: int) -> list[dict]:
    result: list[dict] = []
    for h, d, k in product((0, 1), repeat=3):
        subset = frame[(frame["H"] == h) & (frame["D"] == d) & (frame[f"K_{radius}"] == k)]
        result.append({
            "H": h, "D": d, "K": k,
            "raw_women": int(len(subset)),
            "pnc_incomplete": int((subset["pnc_state"] == "incomplete").sum()),
            "pnc_completed_within_2_days": int((subset["pnc_state"] == "completed").sum()),
            "pnc_unknown_or_dk": int((subset["pnc_state"] == "unknown_or_dk").sum()),
            "independent_clusters": int(subset["point_index"].nunique()),
        })
    return result


def c_k_table(frame: pd.DataFrame, radius: int) -> list[dict]:
    result: list[dict] = []
    for c, k in product((0, 1), repeat=2):
        subset = frame[(frame["C"] == c) & (frame[f"K_{radius}"] == k)]
        result.append({"C": c, "K": k, "raw_women": int(len(subset)), "independent_clusters": int(subset["point_index"].nunique())})
    return result


def audit_country(work: Path, country: str) -> dict:
    points, cluster_to_point = protected_ge(work, country)
    era_months, t2m = load_npz(work / "staging" / "era5_land" / country / "t2m_monthly_values_only.npz", "t2m_kelvin")
    # Assemble the resumable protected CHIRPS vectors in explicit calendar order.
    start, finish = ("1991-01", "2024-03") if country == "nigeria" else ("1991-01", "2021-10")
    chirps_months = [add_month(start, i) for i in range((int(finish[:4]) - int(start[:4])) * 12 + int(finish[-2:]) - int(start[-2:]) + 1)]
    precip_columns: list[np.ndarray] = []
    for value in chirps_months:
        path = work / "external_raw" / "chirps_v3" / country / "monthly_values" / f"chirps_v3_{value.replace('-', '_')}_values_only.npz"
        with np.load(path) as item:
            precip_columns.append(np.asarray(item["precipitation_mm"], dtype="float64"))
    precip = np.column_stack(precip_columns)
    if len(points) != t2m.shape[0] or len(points) != precip.shape[0]:
        raise RuntimeError("Protected point count differs across DHS GE / ERA5-Land / CHIRPS")
    final_needed = add_month(COUNTRIES[country][1], -1)
    if era_months[:1] != ["1991-01"] or era_months[-1] < final_needed or chirps_months[-1] < final_needed:
        raise RuntimeError("Climate period does not cover frozen outcome months")
    heat_monthly, heat_definition = heat(t2m, era_months)
    spi_monthly, spi_definition = spi3(precip, chirps_months)

    nr_file = next((work / "staging" / "official_dhs_members" / country / "nr").glob("*.dta"))
    columns = BASE + ["p3"]
    nr = pd.read_stata(nr_file, columns=columns, convert_categoricals=False)
    _, frame = latest_most_recent_birth(nr)
    early, category = mother_pnc_within_2_days(frame)
    frame = frame.copy()
    frame["outcome_month"] = frame["p3"].map(cmc_to_month)
    expected_start, expected_end = COUNTRIES[country]
    if frame["outcome_month"].min() != expected_start or frame["outcome_month"].max() != expected_end:
        raise RuntimeError("Frozen outcome month coverage differs from prior DHS gate")
    frame["point_index"] = frame["v001"].map(lambda x: cluster_to_point.get(cluster_key(x), np.nan))
    frame["pnc_state"] = np.where(category.isna(), "unknown_or_dk", np.where(early, "completed", "incomplete"))
    usable = frame[frame["point_index"].notna()].copy()
    usable["point_index"] = usable["point_index"].astype(int)
    outcome_months = sorted(usable["outcome_month"].unique())
    event_file = work / "staging" / "ucdp_ged261" / country / "ged261_events_frozen_window.csv"
    events = pd.read_csv(event_file, usecols=["id", "date_start", "date_end", "latitude", "longitude"])
    events["date_start"] = pd.to_datetime(events["date_start"], errors="raise")
    events["date_end"] = pd.to_datetime(events["date_end"], errors="raise")
    if events[["latitude", "longitude"]].isna().any().any():
        raise RuntimeError("UCDP selected event has missing coordinates")
    k_values, supporting_events = k_by_radius(points, events, outcome_months)
    era_index, chirps_index, outcome_index = ({m: i for i, m in enumerate(era_months)}, {m: i for i, m in enumerate(chirps_months)}, {m: i for i, m in enumerate(outcome_months)})
    h_list: list[int] = []; d_list: list[int] = []
    for row in usable.itertuples(index=False):
        prior = [add_month(row.outcome_month, delta) for delta in (-3, -2, -1)]
        if any(month not in era_index or month not in chirps_index for month in prior):
            raise RuntimeError("A frozen prior-three-month window falls outside climate coverage")
        h_list.append(int(heat_monthly[row.point_index, [era_index[m] for m in prior]].any()))
        d_list.append(int(spi_monthly[row.point_index, chirps_index[prior[-1]]] <= -1.0))
    usable["H"] = h_list; usable["D"] = d_list; usable["C"] = (usable["H"].eq(1) & usable["D"].eq(1)).astype(int)
    for radius in RADII:
        usable[f"K_{radius}"] = [int(k_values[radius][row.point_index, outcome_index[row.outcome_month]]) for row in usable.itertuples(index=False)]
    primary = group_table(usable, 50)
    ck_primary = c_k_table(usable, 50)
    structural_zero = [item for item in primary if item["raw_women"] == 0 or item["independent_clusters"] == 0]
    ck_zero = [item for item in ck_primary if item["C"] == 1 and item["K"] == 1 and (item["raw_women"] == 0 or item["independent_clusters"] == 0)]
    return {
        "country": country,
        "outcome_definition": "Latest live/still birth per woman, p19<24; maternal PNC is official DHS-8 m62--m69 timing category 1--3. Newborn PNC is excluded.",
        "outcome_month_coverage": [expected_start, expected_end],
        "prior_window": "Three full calendar months before the outcome month; outcome month excluded.",
        "source_integrity": {
            "era5_land_point_values_sha256": sha256(work / "staging" / "era5_land" / country / "t2m_monthly_values_only.npz"),
            "chirps_monthly_fields": len(chirps_months),
            "ucdp_frozen_event_csv_sha256": sha256(event_file),
            "protected_ge_points": int(len(points)),
            "ge_to_v001_joined_records": int(len(usable)),
            "ge_to_v001_missing_records": int(len(frame) - len(usable)),
        },
        "definitions": {"H": heat_definition, "D": spi_definition, "C": "H=1 and D=1", "K": "At least one UCDP GED 26.1 event overlaps the prior-three-month window and lies within the stated km radius of the protected DHS GE location."},
        "primary_50km": {"h_d_k_groups": primary, "c_k_support": ck_primary, "unique_ucdp_events_supporting_K": supporting_events[50], "structural_zero_h_d_k_groups": structural_zero, "c1_k1_structural_zero": ck_zero},
        "sensitivity_25km": {"status": "PREDECLARED_NOT_EXPANDED_IN_PROOF_PHASE"},
        "sensitivity_100km": {"status": "PREDECLARED_NOT_EXPANDED_IN_PROOF_PHASE"},
        "gate": "FAIL_STOP_STRUCTURAL_SUPPORT" if structural_zero or ck_zero or supporting_events[50] == 0 else "COUNTS_READY_ADEQUACY_THRESHOLD_NOT_PREDECLARED",
        "gate_interpretation": "A structural zero in any predeclared H/D/K cell, a zero C=1/K=1 cell, or no 50-km supporting event requires stopping. No numerical 'sufficient independent clusters/events' minimum was supplied, so this gate does not convert nonzero counts into an adequacy PASS.",
        "guardrail": "Counts only. No association, effect estimate, adjusted model, P value, confidence interval, or public coordinate/cluster identifier is output.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    records = [audit_country(args.work_dir, country) for country in COUNTRIES]
    payload = {"completed_at_utc": datetime.now(timezone.utc).isoformat(), "stage": "external_exposure_support_gate", "records": records, "overall_gate": "FAIL_STOP_STRUCTURAL_SUPPORT" if any(record["gate"] == "FAIL_STOP_STRUCTURAL_SUPPORT" for record in records) else "COUNTS_READY_ADEQUACY_THRESHOLD_NOT_PREDECLARED"}
    target = args.out_dir / "external_exposure_support_manifest.json"
    atomic_json(target, payload)
    atomic_json(args.work_dir / "checkpoints" / "external_exposure_support_gate.complete.json", {"stage": payload["stage"], "completed_at_utc": payload["completed_at_utc"], "artifact": str(target), "overall_gate": payload["overall_gate"]})
    print(json.dumps({"status": "complete", "overall_gate": payload["overall_gate"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
