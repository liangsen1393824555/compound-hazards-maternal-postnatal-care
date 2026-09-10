"""Build the sealed, row-level analysis file from frozen task-02 inputs.

The output is restricted: it contains survey design identifiers and individual
derived records and must remain under work/secure.  This program never writes
coordinates, point climate vectors, event links, or public microdata.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd


SEED = 20260911


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest().upper()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def add_month(value: str, delta: int) -> str:
    year, month = map(int, value.split("-"))
    ordinal = year * 12 + month - 1 + delta
    return f"{ordinal // 12:04d}-{ordinal % 12 + 1:02d}"


def month_start(value: str) -> date:
    year, month = map(int, value.split("-"))
    return date(year, month, 1)


def month_end(value: str) -> date:
    return month_start(add_month(value, 1)) - pd.Timedelta(days=1)


def k_for_windows(points: np.ndarray, events: pd.DataFrame, outcome_months: list[str], deltas: tuple[int, ...], radii: tuple[int, ...]) -> dict[int, np.ndarray]:
    lat1 = np.deg2rad(points[:, 1])[:, None]
    lon1 = np.deg2rad(points[:, 0])[:, None]
    lat2 = np.deg2rad(events["latitude"].to_numpy(dtype="float64"))[None, :]
    lon2 = np.deg2rad(events["longitude"].to_numpy(dtype="float64"))[None, :]
    hav = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    distances = 6371.0088 * 2 * np.arcsin(np.minimum(1, np.sqrt(hav)))
    result = {radius: np.zeros((len(points), len(outcome_months)), dtype=bool) for radius in radii}
    for col, outcome in enumerate(outcome_months):
        first, last = min(deltas), max(deltas)
        start = month_start(add_month(outcome, first))
        end = month_end(add_month(outcome, last))
        active = (events["date_start"].dt.date <= end) & (events["date_end"].dt.date >= start)
        idx = np.flatnonzero(active.to_numpy())
        if len(idx):
            for radius in radii:
                result[radius][:, col] = (distances[:, idx] <= radius).any(axis=1)
    return result


def build_country(upstream: Path, country: str) -> tuple[pd.DataFrame, dict]:
    work = upstream / "work"
    sys.path.insert(0, str(work))
    from external_exposure_support_gate import cmc_to_month, heat, load_npz, protected_ge, spi3
    from internal_variable_denominator_gate import latest_most_recent_birth, mother_pnc_within_2_days

    points, cluster_to_point = protected_ge(work, country)
    era_months, t2m = load_npz(work / "staging" / "era5_land" / country / "t2m_monthly_values_only.npz", "t2m_kelvin")
    climate_start, climate_end = (("1991-01", "2024-03") if country == "nigeria" else ("1991-01", "2021-10"))
    n_months = (int(climate_end[:4]) - int(climate_start[:4])) * 12 + int(climate_end[-2:]) - int(climate_start[-2:]) + 1
    chirps_months = [add_month(climate_start, i) for i in range(n_months)]
    precip_columns = []
    for month in chirps_months:
        p = work / "external_raw" / "chirps_v3" / country / "monthly_values" / f"chirps_v3_{month.replace('-', '_')}_values_only.npz"
        with np.load(p) as item:
            precip_columns.append(np.asarray(item["precipitation_mm"], dtype="float64"))
    precip = np.column_stack(precip_columns)
    heat_monthly, _ = heat(t2m, era_months)
    spi_monthly, _ = spi3(precip, chirps_months)

    nr_file = next((work / "staging" / "official_dhs_members" / country / "nr").glob("*.dta"))
    columns = [
        "v001", "v002", "v003", "v005", "v011", "v021", "v022", "v023", "v024", "v025",
        "v106", "v190", "p3", "p19", "pidx", "pord", "midxp", "m80", "m62", "m63", "m64",
        "m66", "m67", "m68", "m69",
    ]
    nr = pd.read_stata(nr_file, columns=columns, convert_categoricals=False)
    _, frame = latest_most_recent_birth(nr)
    early, category = mother_pnc_within_2_days(frame)
    frame = frame.copy().reset_index(drop=True)
    frame["Y"] = np.where(category.isna(), np.nan, np.where(early, 0, 1))
    frame["outcome_month"] = frame["p3"].map(cmc_to_month)
    frame["outcome_year"] = frame["outcome_month"].str[:4]
    frame["calendar_month"] = frame["outcome_month"].str[-2:]
    frame["point_index"] = frame["v001"].map(lambda value: cluster_to_point.get(str(int(float(value))), np.nan))
    if frame["point_index"].isna().any():
        raise RuntimeError(f"Protected GE link loss in {country}")
    frame["point_index"] = frame["point_index"].astype(int)

    era_index = {m: i for i, m in enumerate(era_months)}
    chirps_index = {m: i for i, m in enumerate(chirps_months)}
    outcomes = sorted(frame["outcome_month"].unique())
    outcome_index = {m: i for i, m in enumerate(outcomes)}
    events_a = pd.read_csv(work / "staging" / "ucdp_ged261" / country / "ged261_events_frozen_window.csv", usecols=["id", "date_start", "date_end", "latitude", "longitude"])
    events_a["date_start"] = pd.to_datetime(events_a["date_start"], errors="raise")
    events_a["date_end"] = pd.to_datetime(events_a["date_end"], errors="raise")
    k_a = k_for_windows(points, events_a, outcomes, (-3, -2, -1), (25, 50, 100))

    h_a: list[int] = []
    d_a: list[int] = []
    ka25: list[int] = []
    ka50: list[int] = []
    ka100: list[int] = []
    for row in frame.itertuples(index=False):
        prior_a = [add_month(row.outcome_month, d) for d in (-3, -2, -1)]
        h_a.append(int(heat_monthly[row.point_index, [era_index[m] for m in prior_a]].any()))
        d_a.append(int(spi_monthly[row.point_index, chirps_index[prior_a[-1]]] <= -1.0))
        oi = outcome_index[row.outcome_month]
        ka25.append(int(k_a[25][row.point_index, oi]))
        ka50.append(int(k_a[50][row.point_index, oi]))
        ka100.append(int(k_a[100][row.point_index, oi]))

    frame["H_A"] = h_a
    frame["D_A"] = d_a
    frame["C_A"] = ((np.asarray(h_a) == 1) & (np.asarray(d_a) == 1)).astype(int)
    frame["K_A25"] = ka25
    frame["K_A50"] = ka50
    frame["K_A100"] = ka100
    for name, ccol, kcol in (
        ("CK_A25", "C_A", "K_A25"), ("CK_A50", "C_A", "K_A50"),
        ("CK_A100", "C_A", "K_A100"),
    ):
        frame[name] = frame[ccol].astype(str) + frame[kcol].astype(str)

    frame["weight"] = frame["v005"] / 1_000_000.0
    frame["age_outcome"] = (frame["p3"] - frame["v011"]) / 12.0
    frame["age_group"] = pd.cut(frame["age_outcome"], bins=[15, 20, 25, 30, 35, 40, 45, 50], right=False, labels=["15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49"])
    frame["pord_group"] = pd.cut(frame["pord"], bins=[0, 1, 4, np.inf], labels=["1", "2-4", "5+"], include_lowest=False)
    frame["country"] = country
    frame["anon_row"] = np.arange(1, len(frame) + 1)

    kept = [
        "country", "anon_row", "Y", "weight", "v005", "v021", "v022", "v023", "v024", "v025",
        "v106", "v190", "p3", "p19", "pord", "m80", "outcome_month", "outcome_year", "calendar_month",
        "age_group", "pord_group", "H_A", "D_A", "C_A", "K_A25", "K_A50", "K_A100", "CK_A25",
        "CK_A50", "CK_A100",
    ]
    out = frame[kept].copy()
    known = out["Y"].notna()
    required = ["CK_A50", "v005", "v021", "v022", "v024", "v025", "age_group", "pord_group", "v106"]
    cc = known & out[required].notna().all(axis=1)
    manifest = {
        "country": country,
        "source_nr_sha256": sha256(nr_file),
        "eligible_n": int(len(out)),
        "outcome_known_n": int(known.sum()),
        "complete_case_n": int(cc.sum()),
        "unknown_n": int((~known).sum()),
    }
    return out, manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--upstream", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()
    np.random.seed(SEED)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    manifests = []
    for country in ("nigeria", "burkina_faso"):
        frame, item = build_country(args.upstream.resolve(), country)
        target = args.out_dir / f"restricted_{country}.csv"
        frame.to_csv(target, index=False, lineterminator="\n", float_format="%.16g")
        item["restricted_csv_sha256"] = sha256(target)
        manifests.append(item)
    scheme_b = {
        "status": "STOP_ITEM_METHOD_UNDERSPECIFIED",
        "reason": "The frozen SAP specifies months -12 through -4 but does not mechanically specify which SPI-3 endpoint rule defines D across that nine-month window; alternative readings change C classification.",
        "model_fitted": False,
    }
    displacement = {
        "status": "STOP_ITEM_INPUT_NOT_MECHANICALLY_AVAILABLE",
        "reason": "Frozen CHIRPS inputs contain values only at original protected GE points, not gridded fields that can be resampled at 200 perturbed locations; SAP requires H/D/K to be recomputed each time.",
        "approximation_used": False,
    }
    write_json(args.out_dir / "derivation_manifest.json", {"seed": SEED, "countries": manifests, "scheme_b_sensitivity": scheme_b, "displacement_sensitivity": displacement})


if __name__ == "__main__":
    main()
