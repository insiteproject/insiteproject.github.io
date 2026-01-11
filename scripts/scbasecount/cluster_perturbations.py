#!/usr/bin/env python3

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path


CONTROL_PATTERNS = [
    r"\bnone\b",
    r"\bno treatment\b",
    r"\buntreated\b",
    r"\bmock\b",
    r"\bvehicle\b",
    r"\bcontrol\b",
    r"\bbaseline\b",
    r"\bnormal\b",
    r"\bhealthy\b",
    r"\bnot specified\b",
    r"\bunsure\b",
    r"\bn/?a\b",
]

VACCINE_PATTERNS = [
    r"\bvaccin",
    r"\bbnt162b2\b",
    r"\bpfizer\b",
    r"\bmoderna\b",
    r"\bjohnson\b",
    r"\bjnj\b",
    r"\bmrna\b",
    r"\bbooster\b",
]

INFECTION_PATTERNS = [
    r"\binfect",
    r"\bcovid\b",
    r"\bsars\b",
    r"\bhiv\b",
    r"\bhbv\b",
    r"\bhcv\b",
    r"\binfluenza\b",
    r"\bflu\b",
    r"\btuberculosis\b",
    r"\btb\b",
    r"\bmycobacter",
    r"\bcmv\b",
    r"\bebv\b",
    r"\brsv\b",
    r"\bhepatitis\b",
    r"\bmalaria\b",
]

AUTOIMMUNE_PATTERNS = [
    r"\bautoimmune\b",
    r"\blupus\b",
    r"\bsle\b",
    r"\brheumatoid\b",
    r"\bcrohn",
    r"\bulcerative colitis\b",
    r"\bpsoriasis\b",
    r"\bceliac\b",
    r"\btype 1 diabetes\b",
    r"\bms\b",
    r"\bmultiple sclerosis\b",
]

CANCER_PATTERNS = [
    r"\bcancer\b",
    r"\btumor\b",
    r"\bmelanoma\b",
    r"\bleukemia\b",
    r"\blymphoma\b",
    r"\bmyeloma\b",
    r"\bcarcinoma\b",
    r"\bmalignant\b",
]

GENETIC_PATTERNS = [
    r"\bcrispr\b",
    r"\bknockout\b",
    r"\bknockdown\b",
    r"\brnai\b",
    r"\bsirna\b",
    r"\bshrna\b",
    r"\bko\b",
    r"\boverexpress",
    r"\btransduce",
    r"\btransfected\b",
    r"\bengineer",
    r"\bcar-?t\b",
    r"\btcr-?t\b",
]

DRUG_PATTERNS = [
    r"\bdrug\b",
    r"\btreat",
    r"\btherapy\b",
    r"\binhibitor\b",
    r"\bblocker\b",
    r"\bagonist\b",
    r"\bantagonist\b",
    r"\bdexamethasone\b",
    r"\bimatinib\b",
    r"\bcheckpoint\b",
    r"\banti-",
]

CYTOKINE_PATTERNS = [
    r"\binterleukin\b",
    r"\bil-?2\b",
    r"\bil-?6\b",
    r"\bil-?7\b",
    r"\bil-?15\b",
    r"\bifn\b",
    r"\binterferon\b",
    r"\btnf\b",
    r"\btgf\b",
    r"\bcytokine\b",
]

STIMULATION_PATTERNS = [
    r"\bactivation\b",
    r"\bstimulat",
    r"\banti-cd3\b",
    r"\banti-cd28\b",
    r"\bpma\b",
    r"\bionomycin\b",
]

DIET_PATTERNS = [
    r"\bdiet\b",
    r"\bglucose\b",
    r"\bfasting\b",
    r"\bobese\b",
    r"\bmetabolic\b",
]

ENV_PATTERNS = [
    r"\bhypoxia\b",
    r"\boxygen\b",
    r"\bradiation\b",
    r"\bheat shock\b",
]


def normalize(text):
    if text is None:
        return ""
    return re.sub(r"\s+", " ", text.strip().lower())


def match_any(patterns, text):
    return any(re.search(pat, text) for pat in patterns)


def classify_perturbation(raw_text):
    text = normalize(raw_text)
    if not text:
        return "control"

    if match_any(CONTROL_PATTERNS, text):
        return "control"
    if match_any(VACCINE_PATTERNS, text):
        return "vaccination"
    if match_any(INFECTION_PATTERNS, text):
        return "infection"
    if match_any(AUTOIMMUNE_PATTERNS, text):
        return "autoimmune"
    if match_any(CANCER_PATTERNS, text):
        return "cancer"
    if match_any(GENETIC_PATTERNS, text):
        return "genetic"
    if match_any(CYTOKINE_PATTERNS, text):
        return "cytokine"
    if match_any(STIMULATION_PATTERNS, text):
        return "stimulation"
    if match_any(DRUG_PATTERNS, text):
        return "drug"
    if match_any(DIET_PATTERNS, text):
        return "diet_metabolic"
    if match_any(ENV_PATTERNS, text):
        return "environment"
    return "other"


def read_summary_csv(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            cleaned = {k: v for k, v in row.items() if k is not None}
            rows.append(cleaned)
    return rows


def write_summary_csv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def update_datasets_json(datasets_path, output_path, summary_map):
    with open(datasets_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    if isinstance(data, list) and data and isinstance(data[0], dict) and "data_objects" in data[0]:
        repos = data
    else:
        repos = [{"data_objects": data}]

    missing = []
    updated = 0

    for repo in repos:
        for obj in repo.get("data_objects", []):
            source_file = None
            extra = obj.get("extra_info") or {}
            if "source_file" in extra:
                source_file = extra["source_file"]
            if not source_file:
                h5ad = obj.get("anndata_object", {}).get("h5ad_path")
                if h5ad:
                    source_file = Path(h5ad).name.replace(".h5ad", ".rds")

            if not source_file or source_file not in summary_map:
                missing.append(source_file or "unknown")
                continue

            meta = summary_map[source_file]
            obj["perturbation_type"] = meta["perturbation_group"]
            obj.setdefault("perturbation_extra_info", {})
            obj["perturbation_extra_info"]["raw"] = meta["perturbation_raw"]
            obj["perturbation_extra_info"]["disease"] = meta.get("disease") or ""
            updated += 1

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(data, f)

    return updated, missing


def main():
    parser = argparse.ArgumentParser(description="Cluster scBaseCount perturbation strings into groups.")
    parser.add_argument("--summary", default="../scbasecount - data/tcells/perturbation_summary.csv",
                        help="CSV produced by extract_perturbations.R")
    parser.add_argument("--out", default="../scbasecount - data/tcells/perturbation_grouped.csv",
                        help="Output CSV with perturbation_group")
    parser.add_argument("--datasets-json", default=None,
                        help="Path to datasets.json to update")
    parser.add_argument("--datasets-out", default=None,
                        help="Output path for updated datasets.json")
    args = parser.parse_args()

    rows = read_summary_csv(args.summary)
    for row in rows:
        raw = row.get("perturbation_raw", "")
        row["perturbation_group"] = classify_perturbation(raw)

    fieldnames = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    if "perturbation_group" not in fieldnames:
        fieldnames.append("perturbation_group")
    write_summary_csv(args.out, rows, fieldnames)

    counts = Counter(r["perturbation_group"] for r in rows)
    print("Perturbation group counts:")
    for key, value in counts.most_common():
        print(f"{key}: {value}")

    if args.datasets_json:
        datasets_out = args.datasets_out or args.datasets_json
        summary_map = {r["file"]: r for r in rows if r.get("file")}
        updated, missing = update_datasets_json(args.datasets_json, datasets_out, summary_map)
        print(f"Updated {updated} data objects.")
        if missing:
            print(f"Missing summary for {len(set(missing))} files (example: {missing[0]})")


if __name__ == "__main__":
    main()
