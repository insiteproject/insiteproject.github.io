#!/usr/bin/env python3

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Dict, List

import numpy as np
import scanpy as sc
import scipy.sparse as sp


def map_tissue_to_location(tissue: str) -> str:
    tissue_lower = tissue.lower() if tissue else "unknown"
    tissue_mapping = {
        "blood": "blood",
        "peripheral blood": "blood",
        "pbmc": "blood",
        "lung": "lung",
        "liver": "liver",
        "brain": "brain",
        "skin": "skin",
        "skin of body": "skin",
        "small intestine": "intestine",
        "large intestine": "intestine",
        "colon": "colon",
        "spleen": "spleen",
        "lymph node": "lymph-node",
        "bone marrow": "bone-marrow",
        "kidney": "kidney",
        "heart": "heart",
        "pancreas": "pancreas",
        "thymus": "thymus",
        "tonsil": "tonsil",
        "tumor": "tumor",
    }
    for key, value in tissue_mapping.items():
        if key in tissue_lower:
            return value
    return tissue_lower.replace(" ", "-")


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


def map_final_type_to_subtype(final_type: str) -> str:
    type_mapping = {
        "CD4": "CD4",
        "CD8": "CD8",
        "Treg": "Treg",
        "T_unspecified": "unknown",
    }
    return type_mapping.get(final_type, "unknown")


def infer_subtype_from_obs(adata: sc.AnnData, column: str) -> str:
    if column not in adata.obs.columns:
        return "unknown"
    counts = adata.obs[column].value_counts()
    for label in counts.index:
        mapped = map_label_to_subtype(label)
        if mapped != "unknown":
            return mapped
    return "unknown"


def infer_donor_type(cell_line: str, tissue: str) -> str:
    if not cell_line:
        return "diseased"
    cell_line_lower = str(cell_line).lower()
    if any(term in cell_line_lower for term in ["healthy", "normal", "control"]):
        return "healthy"
    return "diseased"


def compute_mean_expression_profile(adata: sc.AnnData, hvg_list: List[str]) -> List[List[float]]:
    if adata.n_obs == 0:
        return [[0.0, 0.0] for _ in hvg_list]

    var_index = {g: i for i, g in enumerate(adata.var_names)}
    present_genes = [g for g in hvg_list if g in var_index]

    if not present_genes:
        return [[0.0, 0.0] for _ in hvg_list]

    indices = [var_index[g] for g in present_genes]
    X = adata.X[:, indices]

    if sp.issparse(X):
        means = np.array(X.mean(axis=0)).ravel()
        second_moment = np.array(X.power(2).mean(axis=0)).ravel()
        stds = np.sqrt(np.maximum(second_moment - np.power(means, 2), 0))
    else:
        means = np.mean(X, axis=0)
        stds = np.std(X, axis=0)

    stats = {g: (float(m), float(s)) for g, m, s in zip(present_genes, means, stds)}
    return [[stats[g][0], stats[g][1]] if g in stats else [0.0, 0.0] for g in hvg_list]


def subsample_umap(adata: sc.AnnData, n_points: int = 1000) -> List[List[float]]:
    if adata.n_obs == 0:
        return []

    if adata.n_obs > n_points:
        indices = np.random.choice(adata.n_obs, n_points, replace=False)
        adata_sub = adata[indices].copy()
    else:
        adata_sub = adata.copy()

    try:
        sc.pp.normalize_total(adata_sub, target_sum=1e4)
        sc.pp.log1p(adata_sub)

        if adata_sub.n_vars > 2000:
            sc.pp.highly_variable_genes(adata_sub, n_top_genes=2000, flavor="seurat_v3", span=0.3)
            adata_sub = adata_sub[:, adata_sub.var.highly_variable].copy()

        sc.pp.scale(adata_sub, max_value=10)
        sc.tl.pca(adata_sub, n_comps=min(50, adata_sub.n_vars - 1, adata_sub.n_obs - 1))
        sc.pp.neighbors(adata_sub, n_neighbors=min(15, adata_sub.n_obs - 1))
        sc.tl.umap(adata_sub)

        umap_coords = adata_sub.obsm["X_umap"]
        return [[float(x), float(y)] for x, y in umap_coords]
    except Exception as e:
        print(f"  Warning: UMAP computation failed: {e}")
        return []


def build_datasets(h5ads_dir: Path, hvg_list: List[str], n_umap_points: int) -> Dict:
    data_objects = []
    skipped_mouse = 0
    h5ad_files = sorted(h5ads_dir.glob("*.h5ad"))

    for i, h5ad_file in enumerate(h5ad_files, start=1):
        print(f"[{i}/{len(h5ad_files)}] {h5ad_file.name}")
        adata = sc.read_h5ad(h5ad_file)

        organism = ""
        if "organism" in adata.obs.columns:
            organism = str(adata.obs["organism"].iloc[0])
        if organism and organism != "Homo sapiens":
            skipped_mouse += 1
            continue

        tissue = adata.obs["tissue"].iloc[0] if "tissue" in adata.obs.columns else "unknown"
        cell_line = adata.obs["cell_line"].iloc[0] if "cell_line" in adata.obs.columns else None

        dominant_type = "unknown"
        type_distribution = {}
        if "final_type" in adata.obs.columns:
            type_counts = adata.obs["final_type"].value_counts()
            type_distribution = type_counts.to_dict()
            if len(type_counts) > 0:
                dominant_type = type_counts.index[0]

        location = map_tissue_to_location(tissue)
        t_cell_subtype = map_final_type_to_subtype(dominant_type)
        if t_cell_subtype == "unknown":
            t_cell_subtype = infer_subtype_from_obs(adata, "SingleR_label")
        if t_cell_subtype == "unknown":
            t_cell_subtype = infer_subtype_from_obs(adata, "tc_subtype")
        donor_type = infer_donor_type(cell_line, tissue)

        mean_expr = compute_mean_expression_profile(adata, hvg_list)
        umap_data = subsample_umap(adata, n_points=min(n_umap_points, adata.n_obs))

        data_objects.append({
            "perturbation_type": "control",
            "perturbation_extra_info": {},
            "t_cell_subtype": t_cell_subtype,
            "t_cell_subtype_extra_info": {"distribution": type_distribution},
            "donor_type": donor_type,
            "donor_extra_info": {},
            "time_point": "unknown",
            "location": location,
            "cell_count": adata.n_obs,
            "extra_info": {
                "source_file": f"{h5ad_file.stem}.rds",
                "tissue_original": tissue
            },
            "anndata_object": {
                "h5ad_path": f"h5ads_unified_genes/{h5ad_file.name}",
                "var_columns": ["gene_id"],
                "obs_columns": list(adata.obs.columns),
                "n_obs": adata.n_obs,
                "n_vars": adata.n_vars
            },
            "mean_expression_profile": mean_expr,
            "subsampled_umap": umap_data
        })

    if skipped_mouse:
        print(f"Skipped {skipped_mouse} non-human datasets based on organism metadata.")

    return [{
        "id_on_repo": "scbasecount-tcells",
        "repo_name": "scbasecount",
        "record_name": "T-cell Corpus",
        "url": "https://github.com/scbasecount",
        "description": "T-cell single-cell RNA-seq data from the scBaseCount database. Contains CD4+, CD8+, and Treg T-cells from various tissues including blood, lung, liver, and more.",
        "data_objects": data_objects
    }]


def main():
    parser = argparse.ArgumentParser(description="Build datasets.json from existing h5ad files")
    parser.add_argument("--h5ads-dir", required=True, help="Directory containing h5ads_unified_genes")
    parser.add_argument("--hvg-list", required=True, help="Path to highly_variable_genes_list.json")
    parser.add_argument("--output", required=True, help="Output datasets.json path")
    parser.add_argument("--umap-points", type=int, default=1000, help="Max UMAP points per dataset")
    args = parser.parse_args()

    h5ads_dir = Path(args.h5ads_dir)
    hvg_path = Path(args.hvg_list)
    output_path = Path(args.output)

    if not h5ads_dir.exists():
        raise FileNotFoundError(f"h5ads dir not found: {h5ads_dir}")
    if not hvg_path.exists():
        raise FileNotFoundError(f"HVG list not found: {hvg_path}")

    with hvg_path.open() as f:
        hvg_list = json.load(f)

    np.random.seed(42)
    datasets = build_datasets(h5ads_dir, hvg_list, args.umap_points)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as f:
        json.dump(datasets, f)

    print(f"Saved datasets.json to {output_path}")


if __name__ == "__main__":
    main()
