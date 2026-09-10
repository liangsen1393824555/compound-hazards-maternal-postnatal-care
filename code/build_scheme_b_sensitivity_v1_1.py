"""Derive only SAP-v1.1 Scheme-B exposures in a restricted clean run."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import zipfile
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd


SEED = 20260911
COUNTRY_NAME = {"nigeria": "Nigeria", "burkina_faso": "Burkina Faso"}
ENVELOPES = {
    "nigeria": ("2020-12-01", "2023-12-31"),
    "burkina_faso": ("2018-08-01", "2021-07-31"),
}


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
    y, m = map(int, value.split("-"))
    return date(y, m, 1)


def month_end(value: str) -> date:
    return month_start(add_month(value, 1)) - pd.Timedelta(days=1)


def filtered_scheme_b_events(upstream: Path, country: str, out_path: Path) -> pd.DataFrame:
    archive = upstream / "work" / "external_raw" / "ucdp" / "ged261-csv.zip"
    if sha256(archive) != "8C941D84954E555EE2E54F40FA04D9203BF1E2F962203D0A9930966C4947C667":
        raise RuntimeError("Frozen UCDP global ZIP hash mismatch")
    start, end = ENVELOPES[country]
    selected: list[dict[str, str]] = []
    fields: list[str] | None = None
    with zipfile.ZipFile(archive) as zf:
        if zf.testzip() is not None:
            raise RuntimeError("UCDP ZIP CRC failure")
        member = zf.infolist()[0]
        text = (line.decode("utf-8-sig") for line in zf.open(member))
        reader = csv.DictReader(text)
        fields = reader.fieldnames
        for row in reader:
            if row["country"] == COUNTRY_NAME[country] and row["date_start"][:10] <= end and row["date_end"][:10] >= start:
                selected.append(row)
    if not fields:
        raise RuntimeError("UCDP CSV header missing")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        writer.writeheader(); writer.writerows(selected)
    result = pd.DataFrame([{k: row[k] for k in ("id", "date_start", "date_end", "latitude", "longitude")} for row in selected])
    result["date_start"] = pd.to_datetime(result["date_start"], errors="raise")
    result["date_end"] = pd.to_datetime(result["date_end"], errors="raise")
    result["latitude"] = pd.to_numeric(result["latitude"], errors="raise")
    result["longitude"] = pd.to_numeric(result["longitude"], errors="raise")
    return result


def k_scheme_b(points: np.ndarray, events: pd.DataFrame, outcome_months: list[str]) -> np.ndarray:
    lat1 = np.deg2rad(points[:, 1])[:, None]; lon1 = np.deg2rad(points[:, 0])[:, None]
    lat2 = np.deg2rad(events["latitude"].to_numpy(dtype="float64"))[None, :]
    lon2 = np.deg2rad(events["longitude"].to_numpy(dtype="float64"))[None, :]
    hav = np.sin((lat2-lat1)/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin((lon2-lon1)/2)**2
    distance = 6371.0088 * 2 * np.arcsin(np.minimum(1, np.sqrt(hav)))
    output = np.zeros((len(points), len(outcome_months)), dtype=bool)
    for col, outcome in enumerate(outcome_months):
        start, end = month_start(add_month(outcome, -12)), month_end(add_month(outcome, -4))
        active = (events["date_start"].dt.date <= end) & (events["date_end"].dt.date >= start)
        idx = np.flatnonzero(active.to_numpy())
        if len(idx): output[:, col] = (distance[:, idx] <= 50).any(axis=1)
    return output


def build_country(upstream: Path, country: str, out_dir: Path) -> dict:
    work = upstream / "work"
    sys.path.insert(0, str(work))
    from external_exposure_support_gate import cmc_to_month, heat, load_npz, protected_ge, spi3
    from internal_variable_denominator_gate import latest_most_recent_birth, mother_pnc_within_2_days

    points, cluster_to_point = protected_ge(work, country)
    era_months, t2m = load_npz(work / "staging" / "era5_land" / country / "t2m_monthly_values_only.npz", "t2m_kelvin")
    start, finish = (("1991-01", "2024-03") if country == "nigeria" else ("1991-01", "2021-10"))
    n_months = (int(finish[:4])-int(start[:4]))*12 + int(finish[-2:])-int(start[-2:]) + 1
    chirps_months = [add_month(start, i) for i in range(n_months)]
    precip_cols = []
    for month in chirps_months:
        p = work / "external_raw" / "chirps_v3" / country / "monthly_values" / f"chirps_v3_{month.replace('-', '_')}_values_only.npz"
        with np.load(p) as item: precip_cols.append(np.asarray(item["precipitation_mm"], dtype="float64"))
    heat_monthly, _ = heat(t2m, era_months)
    spi_monthly, _ = spi3(np.column_stack(precip_cols), chirps_months)
    era_index = {m:i for i,m in enumerate(era_months)}; chirps_index = {m:i for i,m in enumerate(chirps_months)}

    nr_file = next((work / "staging" / "official_dhs_members" / country / "nr").glob("*.dta"))
    cols = ["v001","v002","v003","v005","v011","v021","v022","v023","v024","v025","v106","v190","p3","p19","pidx","pord","midxp","m80","m62","m63","m64","m66","m67","m68","m69"]
    nr = pd.read_stata(nr_file, columns=cols, convert_categoricals=False)
    _, frame = latest_most_recent_birth(nr); early, category = mother_pnc_within_2_days(frame)
    frame = frame.copy().reset_index(drop=True)
    frame["Y"] = np.where(category.isna(), np.nan, np.where(early, 0, 1))
    frame["outcome_month"] = frame["p3"].map(cmc_to_month)
    frame["outcome_year"] = frame["outcome_month"].str[:4]; frame["calendar_month"] = frame["outcome_month"].str[-2:]
    frame["point_index"] = frame["v001"].map(lambda x: cluster_to_point.get(str(int(float(x))), np.nan))
    if frame["point_index"].isna().any(): raise RuntimeError(f"GE link loss: {country}")
    frame["point_index"] = frame["point_index"].astype(int)
    outcomes = sorted(frame["outcome_month"].unique()); outcome_index = {m:i for i,m in enumerate(outcomes)}
    event_path = out_dir / f"ucdp_scheme_b_{country}.csv"
    events = filtered_scheme_b_events(upstream, country, event_path)
    kval = k_scheme_b(points, events, outcomes)
    hb=[]; db=[]; kb=[]
    for row in frame.itertuples(index=False):
        w_h = [add_month(row.outcome_month, d) for d in range(-12, -3)]
        w_d = [add_month(row.outcome_month, d) for d in range(-10, -3)]
        hb.append(int(heat_monthly[row.point_index, [era_index[m] for m in w_h]].any()))
        db.append(int((spi_monthly[row.point_index, [chirps_index[m] for m in w_d]] <= -1.0).any()))
        kb.append(int(kval[row.point_index, outcome_index[row.outcome_month]]))
    frame["H_B"] = hb; frame["D_B"] = db
    frame["C_B"] = ((np.asarray(hb)==1) & (np.asarray(db)==1)).astype(int); frame["K_B"] = kb
    frame["CK_B"] = frame["C_B"].astype(str) + frame["K_B"].astype(str)
    frame["weight"] = frame["v005"] / 1_000_000
    age = (frame["p3"]-frame["v011"])/12
    frame["age_group"] = pd.cut(age, [15,20,25,30,35,40,45,50], right=False, labels=["15-19","20-24","25-29","30-34","35-39","40-44","45-49"])
    frame["pord_group"] = pd.cut(frame["pord"], [0,1,4,np.inf], labels=["1","2-4","5+"], include_lowest=False)
    frame["country"] = country; frame["anon_row"] = np.arange(1, len(frame)+1)
    kept = ["country","anon_row","Y","weight","v021","v022","v023","v024","v025","v106","v190","m80","outcome_month","outcome_year","calendar_month","age_group","pord_group","H_B","D_B","C_B","K_B","CK_B"]
    target = out_dir / f"restricted_scheme_b_{country}.csv"
    frame[kept].to_csv(target, index=False, lineterminator="\n", float_format="%.16g")
    return {
        "country": country, "relative_months_H_K": list(range(-12,-3)), "spi3_endpoints_D": list(range(-10,-3)),
        "eligible_n": len(frame), "outcome_known_n": int(frame["Y"].notna().sum()),
        "event_file": event_path.name, "event_rows": len(events), "event_sha256": sha256(event_path),
        "restricted_file": target.name, "restricted_sha256": sha256(target),
    }


def main() -> None:
    p=argparse.ArgumentParser(); p.add_argument("--upstream",type=Path,required=True); p.add_argument("--out-dir",type=Path,required=True)
    args=p.parse_args(); args.out_dir.mkdir(parents=True,exist_ok=True); np.random.seed(SEED)
    records=[build_country(args.upstream.resolve(), c, args.out_dir) for c in ("nigeria","burkina_faso")]
    write_json(args.out_dir/"scheme_b_derivation_manifest.json", {"sap_amendment":"v1.1","seed":SEED,"ucdp_global_sha256":"8C941D84954E555EE2E54F40FA04D9203BF1E2F962203D0A9930966C4947C667","records":records,"displacement_200":"WITHDRAWN_WITHOUT_REPLACEMENT"})


if __name__ == "__main__": main()
