#!/usr/bin/env python3

import argparse
import json
from collections import Counter
from pathlib import Path

import scanpy as sc


def map_label_to_subtype(label: str) -> str:
    if not label:
        return "unknown"
    label_lower = str(label).lower()
    if "treg" in label_lower or "regulatory t" in label_lower:
        return "Treg"
    if "cd4" in label_lower:
        return "CD4"
    if "cd8" in label_lower:
        return "CD8"
    return "unknown"


def infer_subtype_from_obs(adata: sc.AnnData, column: str) -> str:
    if column not in adata.obs.columns:
        return "unknown"
    counts = adata.obs[column].value_counts()
    for label in counts.index:
        mapped = map_label_to_subtype(label)
        if mapped != "unknown":
            return mapped
    return "unknown"


def collapse_donor_type(donor_type: str) -> str:
    if donor_type and donor_type.lower() == "healthy":
        return "healthy"
    return "diseased"


def update_datasets(datasets_path: Path, output_path: Path, h5ads_dir: Path, hvg_limit: int):
    with datasets_path.open("r", encoding="utf-8") as f:
        data = json.load(f)

    if isinstance(data, list) and data and isinstance(data[0], dict) and "data_objects" in data[0]:
        repos = data
    else:
        repos = [{"data_objects": data}]

    updated_subtypes = 0
    skipped_missing = 0

    for repo in repos:
        for obj in repo.get("data_objects", []):
            profile = obj.get("mean_expression_profile")
            if isinstance(profile, list) and len(profile) > hvg_limit:
                obj["mean_expression_profile"] = profile[:hvg_limit]

            obj["donor_type"] = collapse_donor_type(obj.get("donor_type"))

            if obj.get("t_cell_subtype") in (None, "", "unknown"):
                h5ad_rel = obj.get("anndata_object", {}).get("h5ad_path")
                if not h5ad_rel:
                    skipped_missing += 1
                    continue
                h5ad_path = h5ads_dir / Path(h5ad_rel).name
                if not h5ad_path.exists():
                    skipped_missing += 1
                    continue

                adata = sc.read_h5ad(h5ad_path, backed="r")
                try:
                    subtype = infer_subtype_from_obs(adata, "SingleR_label")
                    if subtype == "unknown":
                        subtype = infer_subtype_from_obs(adata, "tc_subtype")
                finally:
                    if hasattr(adata, "file") and adata.file is not None:
                        adata.file.close()
                if subtype != "unknown":
                    obj["t_cell_subtype"] = subtype
                    updated_subtypes += 1

    with output_path.open("w", encoding="utf-8") as f:
        json.dump(data, f)

    return updated_subtypes, skipped_missing


def main():
    parser = argparse.ArgumentParser(description="Shrink HVGs and update donor/t-cell labels for website datasets.")
    parser.add_argument("--datasets", required=True, help="Path to datasets.json")
    parser.add_argument("--output", required=True, help="Output datasets.json path")
    parser.add_argument("--h5ads-dir", required=True, help="Directory containing h5ads_unified_genes")
    parser.add_argument("--hvg-in", required=True, help="Path to existing highly_variable_genes_list.json")
    parser.add_argument("--hvg-out", required=True, help="Path for updated highly_variable_genes_list.json")
    parser.add_argument("--hvg-limit", type=int, default=500, help="Number of HVGs to keep")
    args = parser.parse_args()

    hvg_in = Path(args.hvg_in)
    hvg_out = Path(args.hvg_out)
    with hvg_in.open("r", encoding="utf-8") as f:
        hvgs = json.load(f)

    if not isinstance(hvgs, list) or not hvgs:
        raise ValueError("HVG list is empty or invalid.")

    hvgs = hvgs[: args.hvg_limit]
    hvg_out.parent.mkdir(parents=True, exist_ok=True)
    with hvg_out.open("w", encoding="utf-8") as f:
        json.dump(hvgs, f)

    updated, skipped = update_datasets(
        Path(args.datasets),
        Path(args.output),
        Path(args.h5ads_dir),
        args.hvg_limit,
    )

    print(f"Updated t_cell_subtype on {updated} data objects.")
    if skipped:
        print(f"Skipped {skipped} objects due to missing h5ad path.")


if __name__ == "__main__":
    main()
