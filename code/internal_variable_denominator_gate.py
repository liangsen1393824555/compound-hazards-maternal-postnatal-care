"""Reproducible DHS-8 internal variable and denominator gate.

This program is intentionally limited to official DHS recode data already
verified at intake.  It produces aggregate mappings, denominator flows,
missing/DK counts, and official-indicator reproducibility checks.  It neither
opens external exposure data nor estimates any exposure--outcome association.
"""
from __future__ import annotations

import argparse
import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

IDS = ["v001", "v002", "v003"]
BASE = IDS + ["v005", "p19", "pidx", "midxp", "m80", "m14", "m15", "m62", "m63", "m64", "m66", "m67", "m68", "m69"]
M2 = [f"m2{x}" for x in "abcdefghijklmn"]
SERVICE = ["m80", "m13", "m14", "m15", "m62", "m63", "m64", "m66", "m67", "m68", "m69"]
OFFICIAL = {
    "nigeria": {"anc_skilled": 62.5, "anc4": 52.3, "mother_pnc_2d": 42.9},
    "burkina_faso": {"anc_skilled": 98.3, "anc4": 71.8, "mother_pnc_2d": 78.7},
}


def atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
        temporary = Path(handle.name)
    os.replace(temporary, path)


def weighted_n(frame: pd.DataFrame) -> float:
    return float((frame["v005"] / 1_000_000).sum())


def rate(frame: pd.DataFrame, flag: pd.Series) -> float:
    weights = frame["v005"] / 1_000_000
    return 100 * float(weights[flag.fillna(False)].sum()) / float(weights.sum())


def count_weight(frame: pd.DataFrame, flag: pd.Series) -> dict:
    subset = frame[flag.fillna(False)]
    return {"unweighted_n": int(len(subset)), "weighted_n": round(weighted_n(subset), 1)}


def latest_most_recent_birth(nr: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """DHS standard maternal base: p19 < 24, most recent live/still birth per woman."""
    eligible = nr[(nr["p19"] < 24) & nr["m80"].isin([1, 3])].copy()
    eligible = eligible.sort_values(IDS + ["pidx"], kind="stable")
    latest = eligible.drop_duplicates(IDS, keep="first").copy()
    return eligible, latest


def mother_pnc_within_2_days(frame: pd.DataFrame) -> pd.Series:
    """DHS-8 RH_PNC.do maternal timing logic (lines 47-109), categories 1--3."""
    # 999 means "PNC reported, timing unresolved".  The sequence mirrors the
    # official program: before-discharge timing takes precedence, then after-
    # discharge/home timing; no reported check is assigned zero.
    timing = pd.Series(float("nan"), index=frame.index, dtype="float64")
    momcheck = (frame["m62"].eq(1) | frame["m66"].eq(1))
    timing.loc[momcheck] = 999
    before = frame["m64"].between(11, 29, inclusive="both")
    timing.loc[before] = frame.loc[before, "m63"]
    before_immediate = (timing < 1000) & frame["m64"].between(31, 99, inclusive="both")
    timing.loc[before_immediate] = 0
    after = timing.eq(999) & frame["m68"].between(11, 29, inclusive="both")
    timing.loc[after] = frame.loc[after, "m67"]
    after_immediate = timing.eq(999) & (frame["m67"] < 1000) & frame["m68"].between(31, 99, inclusive="both")
    timing.loc[after_immediate] = 0
    timing.loc[~momcheck] = 0
    category = pd.Series(float("nan"), index=frame.index, dtype="float64")
    category.loc[timing.isin(list(range(0, 1)) + list(range(242, 300)) + list(range(306, 900)))] = 0
    category.loc[timing.isin([100, 101, 102, 103])] = 1
    category.loc[timing.isin(list(range(104, 124)) + [200])] = 2
    category.loc[timing.isin(list(range(124, 172)) + [201, 202])] = 3
    category.loc[timing.isin(list(range(172, 198)) + [203, 204, 205, 206])] = 4
    category.loc[timing.isin(list(range(207, 242)) + list(range(301, 306)))] = 5
    return category.isin([1, 2, 3]), category


def compact_state(frame: pd.DataFrame, flag: pd.Series, unknown: pd.Series | None = None) -> dict:
    response = {"yes": count_weight(frame, flag), "no": count_weight(frame, ~flag)}
    if unknown is not None:
        response["unknown_or_dk"] = count_weight(frame, unknown)
        known = ~unknown
        response["no"] = count_weight(frame, known & ~flag)
    return response


def mapping_qc(ir_file: Path, nr: pd.DataFrame) -> dict:
    reader = pd.io.stata.StataReader(ir_file)
    available = list(reader.variable_labels())
    columns = IDS + [f"{variable}_{suffix}" for variable in SERVICE for suffix in range(1, 7) if f"{variable}_{suffix}" in available]
    ir = pd.read_stata(ir_file, columns=columns, convert_categoricals=False)
    pieces = []
    for suffix in range(1, 7):
        selection = IDS + [f"{variable}_{suffix}" for variable in SERVICE if f"{variable}_{suffix}" in ir]
        wide = ir[selection].copy()
        wide["pidx"] = suffix
        wide = wide.rename(columns={f"{variable}_{suffix}": variable for variable in SERVICE if f"{variable}_{suffix}" in wide})
        pieces.append(wide)
    long_ir = pd.concat(pieces, ignore_index=True)
    service_present = [x for x in SERVICE if x in long_ir and x in nr]
    nonempty = long_ir[service_present].notna().any(axis=1)
    left = long_ir[nonempty].copy()
    right = nr[nr["pidx"].between(1, 6) & nr["m80"].notna()][IDS + ["pidx"] + service_present].copy()
    merged = left.merge(right, on=IDS + ["pidx"], how="outer", suffixes=("_ir", "_nr"), indicator=True)
    equality = {}
    complete = merged[merged["_merge"].eq("both")]
    for variable in service_present:
        a, b = complete[f"{variable}_ir"], complete[f"{variable}_nr"]
        equality[variable] = int((a.eq(b) | (a.isna() & b.isna())).sum()) == int(len(complete))
    return {
        "rule": "IR field X_1...X_6 is reshaped to long with pidx=1...6 and linked to NR by v001/v002/v003/pidx; NR is canonical for analysis.",
        "ir_nonempty_service_sections": int(len(left)),
        "nr_outcome_records_pidx_1_to_6": int(len(right)),
        "merged_both": int((merged["_merge"] == "both").sum()),
        "ir_without_nr": int((merged["_merge"] == "left_only").sum()),
        "nr_without_ir": int((merged["_merge"] == "right_only").sum()),
        "all_compared_values_equal": equality,
    }


def audit_country(staging: Path, country: str) -> dict:
    nr_file = next((staging / country / "nr").glob("*.dta"))
    ir_file = next((staging / country / "ir").glob("*.dta"))
    nr = pd.read_stata(nr_file, columns=BASE + M2, convert_categoricals=False)
    all_eligible, latest = latest_most_recent_birth(nr)
    pnc_early, pnc_category = mother_pnc_within_2_days(latest)
    anc4 = latest["m14"].ge(4) & latest["m14"].lt(90)
    anc_dk = latest["m14"].ge(90)
    if country == "nigeria":
        skilled_columns = ["m2a", "m2b"]
    else:
        skilled_columns = ["m2a", "m2b", "m2c", "m2d", "m2e", "m2f"]
    skilled = latest[skilled_columns].eq(1).any(axis=1)
    # DHS place coding: 11/12 are homes; 21--89 are named facilities; 96 is
    # retained as DK/other rather than silently assigned to a facility.
    facility = latest["m15"].between(21, 89, inclusive="both")
    delivery_dk = latest["m15"].ge(90) | latest["m15"].isna()
    pnc_unknown = pnc_category.isna()
    official = OFFICIAL[country]
    reproduced = {
        "anc_skilled": round(rate(latest, skilled), 1),
        "anc4": round(rate(latest, anc4), 1),
        "mother_pnc_2d": round(rate(latest, pnc_early), 1),
    }
    differences = {name: round(reproduced[name] - reference, 1) for name, reference in official.items()}
    return {
        "country": country,
        "scope": "DHS-only internal structural gate; no external data or associations.",
        "mapping_qc": mapping_qc(ir_file, nr),
        "denominator": {
            "definition": "NR p19<24; m80=1 (most recent live birth) or 3 (most recent stillbirth); sort by IDs+pidx and retain one most-recent eligible record per woman.",
            "pre_deduplication": {"unweighted_n": int(len(all_eligible)), "weighted_n": round(weighted_n(all_eligible), 1)},
            "latest_record_per_woman": {"unweighted_n": int(len(latest)), "weighted_n": round(weighted_n(latest), 1), "duplicate_eligible_records_removed": int(len(all_eligible) - len(latest))},
            "outcome_type_pre_deduplication": {str(int(k)): int(v) for k, v in all_eligible["m80"].value_counts().sort_index().items()},
            "outcome_type_latest": {str(int(k)): int(v) for k, v in latest["m80"].value_counts().sort_index().items()},
            "pidx_equals_midxp_all_service_eligible": bool((all_eligible["pidx"] == all_eligible["midxp"]).all()),
        },
        "node_states": {
            "anc_skilled": {"definition": "any country-predeclared qualified-provider m2 item=1", "qualified_provider_fields": skilled_columns, **compact_state(latest, skilled)},
            "anc4": {"definition": "m14=4--89; m14>=90 is DK", **compact_state(latest, anc4, anc_dk)},
            "facility_delivery": {"definition": "m15=21--89 named facility; m15>=90/missing is DK/other", **compact_state(latest, facility, delivery_dk)},
            "mother_pnc_within_2_days": {"definition": "official DHS-8 RH_PNC.do maternal timing categories 1--3; mother variables m62--m69 only", **compact_state(latest, pnc_early, pnc_unknown)},
        },
        "official_weighted_indicator_reproduction": {"reference_percent": official, "reproduced_percent": reproduced, "difference_percentage_points": differences},
        "newborn_separation": "Newborn PNC uses m70--m76 and is deliberately excluded. Mother PNC uses m62--m69 only.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    payload = {
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "gate": "internal_variable_and_denominator",
        "countries": [audit_country(args.work_dir / "staging" / "official_dhs_members", c) for c in ("nigeria", "burkina_faso")],
        "guardrail": "No H/D/K effect, association, model, P value, or external exposure support tabulation is calculated by this program.",
    }
    artifact = args.out_dir / "internal_variable_denominator_gate.json"
    atomic_json(artifact, payload)
    atomic_json(args.work_dir / "checkpoints" / "internal_variable_denominator_gate.complete.json", {"stage": "internal_variable_denominator_gate", "completed_at_utc": payload["completed_at_utc"], "artifact": str(artifact), "resume": "Rerun this program; it replaces only this aggregate checkpoint atomically."})
    print(json.dumps({"status": "complete", "countries": len(payload["countries"])}, ensure_ascii=False))


if __name__ == "__main__":
    main()
