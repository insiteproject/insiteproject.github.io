# INSITE Predict Phase (Observational → Perturbational Prediction)

This phase predicts gene-expression responses to *unseen gene perturbations* using **control-only observational data** from the scBaseCount **T-cell corpus**.

## What’s Implemented Now

- **KO/KD benchmark builder** (offline): `insite/scripts/predict/build_ko_kd_benchmark.py`
- **Method 1 baseline** (offline + website): **Gaussian conditional** predictor trained on control datasets
- **Method 2 (deep)** (offline + website): **DAE denoiser** trained on control pseudobulk means
- **Method 3 (invariant)** (offline + website): **within-BioProject invariant precision** model
- **Website wiring** (Predict tab): `Naive KO` vs `Gaussian (Obs-only)` vs `Deep (DAE)` vs `Invariant`

## Key Data Inputs

- Website summaries (fast, used for UI + model training on HVGs):
  - `insiteproject.github.io/data/datasets.json` (per-sample metadata + mean/std profiles)
  - `insiteproject.github.io/data/highly_variable_genes_list.json` (250 displayed genes; preserves order)
- Full single-cell data (offline only, used to compute means for KO genes not in the 250 HVGs):
  - `tcells_pipeline/output_all_human/h5ads_unified_genes/*.h5ad`

## Benchmark Construction (KO/KD)

Goal: build a clean evaluation set of **(control context → perturbed outcome)** pairs for **gene KO/KD**.

`insite/scripts/predict/build_ko_kd_benchmark.py`:

1. Scans `perturbation_extra_info.raw` to find **gene KO/KD-like** perturbations.
2. Extracts gene symbols (validated against the unified gene universe in the h5ads).
3. Matches a control set by context (t-cell subtype + donor type + tissue + time point), relaxing fields if needed.
4. Computes:
   - Control baseline mean/std for the 250 HVGs (from `datasets.json`, weighted by `cell_count`)
   - KO-gene means for control and perturbed samples (from `.h5ad`, weighted)

Outputs:
- Benchmark JSON (contains arrays): e.g. `insite/scripts/predict/output/ko_kd_benchmark.json`
- Metadata CSV (no arrays): e.g. `insite/scripts/predict/output/ko_kd_benchmark.csv`

Note: the current scBaseCount tcells slice contains **few gene-specific KO/KD samples**; many “CRISPR pool/screen” entries are excluded because they do not specify a target gene in the processed metadata.

## Method 1 (Implemented): Gaussian Conditional Baseline

This method learns gene–gene dependence from controls and predicts intervention effects via a Gaussian conditional.

Training (per cell type, controls only):
- Build a matrix of per-sample mean expression vectors.
- Apply `log1p` to means (stabilizes the count scale).
- Fit a Ledoit–Wolf shrinkage covariance estimator and store the **precision** matrix.

Inference (context-conditioned):
- Compute a context baseline mean `μ` from **control datasets only** (same filters as the UI).
- For KO genes `S`, set `x_S = 0` (in `log1p` space).
- Predict remaining genes with the Gaussian conditional mean using precision submatrices:
  - `E[x_R | x_S] = μ_R - P_RS (P_SS)^-1 (x_S - μ_S)`
- Convert back with `expm1`.

Website artifacts (generated):
- `insiteproject.github.io/data/predict/models/gaussian_conditional_cd4.json`
- `insiteproject.github.io/data/predict/models/gaussian_conditional_cd8.json`
- `insiteproject.github.io/data/predict/models/gaussian_conditional_treg.json`

Website inference code:
- `insiteproject.github.io/js/predict.js`

## Method 2 (Implemented): Denoising Autoencoder (DAE) Baseline

This method fits a small denoising MLP in **log1p space** on control-only pseudobulk means and uses it as a learned gene-dependence denoiser.

Training:
- Input: per-sample `log1p(mean_expression)` vectors (controls only, per cell type)
- Corrupt by random masking + small Gaussian noise
- Train MLP to reconstruct the clean vector (MSE)

KO inference (heuristic):
- Start from the context baseline mean `μ` (controls only)
- Set KO genes to 0 (in log1p space), run one forward pass through the denoiser, clamp KO genes to 0 again
- Convert back with `expm1`

Artifacts:
- `insiteproject.github.io/data/predict/models/dae_denoiser_cd4.json`
- `insiteproject.github.io/data/predict/models/dae_denoiser_cd8.json`
- `insiteproject.github.io/data/predict/models/dae_denoiser_treg.json`

Trainer:
- `insite/scripts/predict/train_dae_baseline.py`

## Method 3 (Implemented): Invariant Precision via Within-BioProject Centering

Goal: estimate gene–gene dependence that is less sensitive to study/batch mean shifts.

Training:
1. Map each sample to a **BioProject** (study-of-origin proxy) via SRA runinfo.
2. For each BioProject (environment), compute the mean log-expression across control datasets.
3. Fit Ledoit–Wolf precision on the **within-environment residuals**.

Inference:
- Same conditional-Gaussian formula as Method 1, but using the invariant precision matrix.

Artifacts:
- `insiteproject.github.io/data/predict/models/invariant_within_bioproject_cd4.json`
- `insiteproject.github.io/data/predict/models/invariant_within_bioproject_cd8.json`
- `insiteproject.github.io/data/predict/models/invariant_within_bioproject_treg.json`

Trainers:
- `insite/scripts/predict/fetch_sra_bioproject_map.py`
- `insite/scripts/predict/train_invariant_precision.py`

## How To Run (Offline)

Build benchmark:

```bash
python3 insite/scripts/predict/build_ko_kd_benchmark.py \
  --datasets insiteproject.github.io/data/datasets.json \
  --genes insiteproject.github.io/data/highly_variable_genes_list.json \
  --h5ads-dir tcells_pipeline/output_all_human/h5ads_unified_genes \
  --out-json insite/scripts/predict/output/ko_kd_benchmark.json \
  --out-csv insite/scripts/predict/output/ko_kd_benchmark.csv
```

Train website baseline models (HVG-only):

```bash
python3 insite/scripts/predict/train_gaussian_baseline.py \
  --datasets insiteproject.github.io/data/datasets.json \
  --genes insiteproject.github.io/data/highly_variable_genes_list.json \
  --out-dir insiteproject.github.io/data/predict/models
```

Train Method 2 (DAE) models (HVG-only):

```bash
python3 insite/scripts/predict/train_dae_baseline.py \
  --datasets insiteproject.github.io/data/datasets.json \
  --genes insiteproject.github.io/data/highly_variable_genes_list.json \
  --out-dir insiteproject.github.io/data/predict/models
```

Fetch BioProject mapping + train Method 3 (Invariant) models (HVG-only):

```bash
python3 insite/scripts/predict/fetch_sra_bioproject_map.py \
  --datasets insiteproject.github.io/data/datasets.json \
  --out insite/scripts/predict/output/sra_bioproject_map.csv

python3 insite/scripts/predict/train_invariant_precision.py \
  --datasets insiteproject.github.io/data/datasets.json \
  --genes insiteproject.github.io/data/highly_variable_genes_list.json \
  --env-map insite/scripts/predict/output/sra_bioproject_map.csv \
  --out-dir insiteproject.github.io/data/predict/models
```

Train + evaluate with KO genes included (offline-only; needs h5ads):

```bash
python3 insite/scripts/predict/train_gaussian_baseline.py \
  --datasets insiteproject.github.io/data/datasets.json \
  --genes insiteproject.github.io/data/highly_variable_genes_list.json \
  --benchmark insite/scripts/predict/output/ko_kd_benchmark.json \
  --h5ads-dir tcells_pipeline/output_all_human/h5ads_unified_genes \
  --out-dir insite/scripts/predict/output/models_with_ko

python3 insite/scripts/predict/evaluate_gaussian_baseline.py \
  --benchmark insite/scripts/predict/output/ko_kd_benchmark.json \
  --model-dir insite/scripts/predict/output/models_with_ko \
  --out insite/scripts/predict/output/eval_gaussian_baseline.json
```

You can also evaluate invariant/DAE variants (with KO genes included in the model gene list):

```bash
python3 insite/scripts/predict/train_invariant_precision.py \
  --datasets insiteproject.github.io/data/datasets.json \
  --genes insiteproject.github.io/data/highly_variable_genes_list.json \
  --env-map insite/scripts/predict/output/sra_bioproject_map.csv \
  --benchmark insite/scripts/predict/output/ko_kd_benchmark.json \
  --h5ads-dir tcells_pipeline/output_all_human/h5ads_unified_genes \
  --out-dir insite/scripts/predict/output/models_invariant_with_ko

python3 insite/scripts/predict/evaluate_precision_model.py \
  --benchmark insite/scripts/predict/output/ko_kd_benchmark.json \
  --model-dir insite/scripts/predict/output/models_invariant_with_ko \
  --model-template 'invariant_within_bioproject_{cell}.json' \
  --out insite/scripts/predict/output/eval_invariant_with_ko.json

python3 insite/scripts/predict/train_dae_baseline.py \
  --datasets insiteproject.github.io/data/datasets.json \
  --genes insiteproject.github.io/data/highly_variable_genes_list.json \
  --benchmark insite/scripts/predict/output/ko_kd_benchmark.json \
  --h5ads-dir tcells_pipeline/output_all_human/h5ads_unified_genes \
  --out-dir insite/scripts/predict/output/models_dae_with_ko

python3 insite/scripts/predict/evaluate_dae_model.py \
  --benchmark insite/scripts/predict/output/ko_kd_benchmark.json \
  --model-dir insite/scripts/predict/output/models_dae_with_ko \
  --out insite/scripts/predict/output/eval_dae_with_ko.json
```

## scPerturb Benchmark + Evaluation (ShifrutMarson2018)

Because scBaseCount has too few clean KO/KD samples, `plans.md` switches Predict evaluation to **scPerturb** (curated Perturb-seq).

Dataset:
- `/fast/datawork/perturbation_data/ShifrutMarson2018.h5ad` (CD8 T cells, CRISPR KO, **Stim** condition has perturbations)

One-command benchmark + training + replicate CV:

```bash
python3 insite/scripts/predict/run_scperturb_shifrut_eval.py \
  --h5ad perturbation_data/ShifrutMarson2018.h5ad \
  --base-genes insiteproject.github.io/data/highly_variable_genes_list.json \
  --out-dir insite/scripts/predict/output/scperturb_shifrut \
  --min-cells 100 \
  --n-samples-per-library 200 \
  --cells-per-sample 200 \
  --dae-damping 0.2
```

Outputs:
- `insite/scripts/predict/output/scperturb_shifrut/benchmark.json` (40 KO groups: 20 genes × 2 replicates, Stim only)
- `insite/scripts/predict/output/scperturb_shifrut/models_all/*.json` (Gaussian / Invariant / DAE models trained on control cells)
- `insite/scripts/predict/output/scperturb_shifrut/eval_*.json` (method eval reports)
- `insite/scripts/predict/output/scperturb_shifrut/cv_test_D{1,2}/` (replicate CV folds)
- `insite/scripts/predict/output/scperturb_shifrut/summary.json` (combined summary)

## Cross-Dataset Evaluation (scBaseCount → scPerturb with Domain Alignment)

Evaluates scBaseCount-trained models on scPerturb T-cell perturbation data, learning an alignment function between controls to handle domain mismatch.

```bash
python3 insite/scripts/predict/eval_cross_dataset_aligned.py \
  --scbasecount-datasets insiteproject.github.io/data/datasets.json \
  --scbasecount-genes insiteproject.github.io/data/highly_variable_genes_list.json \
  --scbasecount-models insiteproject.github.io/data/predict/models \
  --scperturb-h5ad perturbation_data/ShifrutMarson2018.h5ad \
  --out-dir insite/scripts/predict/output/cross_dataset_eval \
  --condition Stim \
  --cell-type CD8
```

Outputs:
- `insite/scripts/predict/output/cross_dataset_eval/cross_dataset_eval.json`

### Key Finding: Observational vs Causal Correlations

The precision-based methods (Gaussian, Invariant) achieve **strong negative correlations** (-0.78) with observed KO effects, but **strong positive correlations** (+0.78) when predictions are negated. This reveals:

> **Observational gene correlations are anti-correlated with causal effects.**

The models successfully identify which genes will change (top-k overlap ~50%), but predict the wrong direction. The knocked-out genes in this dataset are mostly **repressors** - their targets increase (not decrease) when knocked out.

| Method | Pearson | Negated Pearson | Top-k |
|--------|---------|-----------------|-------|
| Gaussian (aligned) | -0.76 | **+0.76** | 0.53 |
| DAE (aligned) | **+0.28** | -0.28 | 0.43 |
| Naive | 0.00 | 0.00 | 0.13 |
