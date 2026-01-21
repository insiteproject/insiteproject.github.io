# Data Flow (scBaseCount tcells -> INSITE)

This document describes how scBaseCount tcells data moves through the pipeline
and into the INSITE website.

## Inputs

- scBaseCount tcells RDS files: `/fast/datawork/scbasecount - data/tcells`
- Selected human-only file list: `/fast/datawork/tcells_pipeline/selected_files_all_human.txt`

## Processing Steps

1) RDS -> h5ad (per dataset)
   - Script: `tcells_pipeline/convert_tcells_to_schema.py`
   - Output: `tcells_pipeline/output_all_human/h5ads_unified_genes/*.h5ad`

2) Dataset metadata + expression profiles
   - Script: `insite/scripts/scbasecount/update_datasets_for_site.py`
   - Actions:
     - Trim HVGs to the configured count (currently 250)
     - Flatten `mean_expression_profile` to reduce JSON size
     - Collapse donor types into `healthy` / `diseased`
     - Fill missing T-cell subtype using `SingleR_label` (fallback: `tc_subtype`)
     - Drop per-dataset `subsampled_umap` (site uses unified UMAP)
   - Outputs:
     - `insite/data/datasets.json`
     - `insite/data/highly_variable_genes_list.json`

3) Unified UMAP (global embedding)
   - Script: `insite/scripts/scbasecount/build_unified_umap.py`
   - Inputs:
     - h5ads
     - `datasets.json` metadata
     - HVG list (250)
   - Steps:
     - Subsample up to 1000 cells per dataset
     - Concatenate all subsamples
     - Normalize + log1p + ComBat batch correction (batch key = `source_file`)
     - PCA -> neighbors -> UMAP
   - Output:
     - `insite/data/unified_umap.json`

## Website Consumption

- Metadata + expression profiles:
  - `insiteproject.github.io/data/datasets.json`
  - `insiteproject.github.io/data/highly_variable_genes_list.json`
- Unified UMAP:
  - `insiteproject.github.io/data/unified_umap.json`
- JS:
  - `insiteproject.github.io/js/explore.js` loads `unified_umap.json` first and
    falls back to per-dataset UMAP if missing.

## Notes

- The unified UMAP is a single shared coordinate space for cross-dataset
  visualization. Per-dataset UMAPs were removed to keep load time fast and to
  avoid mixing incompatible coordinate systems.
