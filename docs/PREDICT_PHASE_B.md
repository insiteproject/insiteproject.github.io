# INSITE Predict Tab - Phase B Implementation

**Date:** 2026-01-20
**Purpose:** Documentation for the Predict tab website deployment (Phase B)

---

## Overview

The Predict tab allows users to simulate gene knockout (KO) effects on T-cell expression profiles using three statistical/ML methods trained on observational (control-only) data from scBaseCount.

**Key Features:**
- Select context filters (T-cell subtype, donor type, tissue)
- Choose 1-5 genes to knock out from the 250 HVG panel
- Compare predictions across 4 methods (Naive, Gaussian, Invariant, DAE)
- View ranked expression changes (delta)
- Export predictions to CSV

---

## Baseline Computation

### Control-Only Filtering

The baseline expression profile is computed from **control-only** datasets:

```javascript
const expressionProfile = await dataService.computeExpressionProfile({
    ...filters,
    perturbation_type: 'control'  // Always filter to controls
});
```

### Weighted Average

For each gene, the baseline is a cell-count-weighted average across matching datasets:

```
mu_counts[gene] = sum(weight_i * mean_i) / sum(weight_i)
```

Where:
- `weight_i` = dataset cell count (`cell_count` or `n_obs`)
- `mean_i` = dataset mean expression from `mean_expression_profile`

### Data Format

`mean_expression_profile` in datasets.json is flattened:
```
[mean_0, std_0, mean_1, std_1, ..., mean_249, std_249]
```

Extract means: `profile[::2]` or `profile[i*2]` for gene i

### Caching

Baseline computations are cached by filter key:
```javascript
const filterKey = JSON.stringify({
    t_cell_subtype, donor_type, time_point, location
});
baselineCache.set(filterKey, { expression, uncertainty });
```

---

## Knockout Simulation

### Input Transformation

1. Convert baseline to log-space: `mu_log = log1p(mu_counts)`
2. Set KO gene values to zero: `mu_log[ko_indices] = 0`

### KO Semantics

KO genes are clamped to zero in log-space (equivalent to zero counts):
- `log1p(0) = 0`
- After expm1: predicted counts = 0

---

## Inference Methods

### Method A: Naive KO

Simply sets KO gene expression to 0, leaves others unchanged.

```javascript
expression[gene] = selectedGenes.has(gene) ? 0 : originalExpression[gene];
```

### Method B: Gaussian Conditional (Precision Model)

Uses a learned precision matrix P from control-only observational data.

**Inference formula:**
Given KO indices S and remaining indices R:

```
x_R = mu_R - P_RS @ solve(P_SS + ridge*I, x_S - mu_S)
```

Where:
- `P_RS` = precision submatrix (remaining × knocked-out)
- `P_SS` = precision submatrix (knocked-out × knocked-out)
- `x_S = 0` (KO genes clamped to zero in log-space)
- `ridge = 1e-6` for numerical stability

**Implementation:**
```javascript
const rhs = fixedIdx.map(i => 0 - mu[i]);
const alpha = solveLinearSystem(Pss, rhs);
for (let i of remainingIndices) {
    pred[i] = mu[i] - dot(P_row_i, alpha);
}
```

### Method C: Invariant Precision (Within-BioProject)

Same inference as Gaussian, but precision matrix learned with within-BioProject centering to remove batch effects.

### Method D: DAE Denoiser (MLP)

Denoising autoencoder trained on control data.

**Forward pass:**
```javascript
x = modelGenes.map(g => log1p(expression[g]));
x[ko_indices] = 0;  // Pre-clamp
y = mlpForward(x, layers);
y[ko_indices] = 0;  // Post-clamp
```

**Activation:** GELU

---

## Output Computation

### Delta Calculation

```javascript
delta[gene] = predicted_counts[gene] - baseline_counts[gene]
```

### Conversion Back to Counts

```javascript
pred_counts = expm1(pred_log)  // Clip to [0, max]
```

---

## Model File Schema

### Precision Models (Gaussian/Invariant)

Location: `data/predict/models/gaussian_conditional_{celltype}.json`

```json
{
    "version": 1,
    "model_type": "gaussian_conditional" | "invariant_within_env_precision",
    "cell_type": "CD4" | "CD8" | "Treg",
    "genes": ["gene1", "gene2", ...],  // Must match HVG order exactly
    "transform": { "name": "log1p" },
    "training": {
        "n_control_datasets": 407,
        "gene_count": 250
    },
    "precision": [[...], [...], ...]  // 250x250 matrix
}
```

### DAE Models

Location: `data/predict/models/dae_denoiser_{celltype}.json`

```json
{
    "version": 1,
    "model_type": "denoising_mlp",
    "cell_type": "CD4" | "CD8" | "Treg",
    "genes": ["gene1", "gene2", ...],
    "transform": { "name": "log1p" },
    "training": {
        "hidden": 128,
        "depth": 2,
        "mask_prob": 0.4,
        "noise_std": 0.07
    },
    "mlp": {
        "layers": [
            { "weight": [[...]], "bias": [...], "activation": "gelu" },
            { "weight": [[...]], "bias": [...], "activation": "linear" }
        ]
    }
}
```

---

## Gene Validation

Models are validated at load time to ensure gene alignment:

```javascript
function validateModelGenes(model, modelName) {
    if (modelGenes.length !== hvgGenes.length) {
        console.error('Gene count mismatch');
        return false;
    }
    // Check first/last genes match HVG list
    for (const i of [0, 1, 2, n-3, n-2, n-1]) {
        if (modelGenes[i] !== hvgGenes[i]) return false;
    }
    return true;
}
```

---

## UI Components

### Gene Search

Filter the 250-gene grid by name substring:
```javascript
function filterGeneGrid(searchTerm) {
    cells.forEach(cell => {
        const matches = geneName.toLowerCase().includes(searchTerm);
        cell.classList.toggle('gene-hidden', !matches);
    });
}
```

### Delta Visualization

1. **Ranked Tables:** Top upregulated and downregulated genes
2. **Heatmap:** Color-coded delta values with diverging scale
3. **Comparison Table:** Side-by-side method comparison

### CSV Export

Exports all 250 genes with:
- Baseline expression
- Predicted expression (all methods)
- Delta values (all methods)
- KO flag

---

## Performance Considerations

1. **Baseline Caching:** Computed once per filter combination
2. **Model Caching:** Loaded once per cell type
3. **Linear Solve:** Gaussian elimination (O(k^3) for k KO genes)
4. **DAE Forward:** Matrix multiplications (O(n*hidden))

Expected latency: < 1 second for typical predictions

---

## Offline Model Training

Models are trained offline using Python scripts in `insite/scripts/predict/`:

```bash
# Train Gaussian model
python train_gaussian_baseline.py \
    --datasets insiteproject.github.io/data/datasets.json \
    --genes insiteproject.github.io/data/highly_variable_genes_list.json \
    --out-dir output/

# Train DAE model
python train_dae_baseline.py --cell-type CD8

# Train Invariant model
python train_invariant_precision.py --env-map output/sra_bioproject_map.csv
```

Copy output JSONs to `insiteproject.github.io/data/predict/models/`

---

## Evaluation Context

See `EVALUATION_FIX_LOG.md` for offline evaluation against scPerturb data:

- **Sign inversion:** Unaligned models show ~-0.51 Pearson with KO effects
- **Alignment helps:** Gaussian/Invariant with domain alignment: +0.34 Pearson
- **DAE struggles:** Alignment doesn't help DAE (-0.49 Pearson)

**Note:** Website uses unaligned models (no scPerturb at runtime). Predictions should be interpreted as "expected downstream effects under control conditions" not "causal perturbation effects."

---

## Known Limitations

1. **250 HVG panel:** Many perturbation targets not in gene set
2. **No alignment at runtime:** Cannot correct for domain shift
3. **Observational learning:** Models learn correlations, not causation
4. **Single-cell averaging:** Pseudobulk loses cell-level heterogeneity

---

## Files Modified (Phase B)

### HTML
- `predict.html` - Added gene search, delta views, export button

### JavaScript
- `js/predict.js` - Added:
  - Gene search/filter
  - Baseline caching
  - Gene validation
  - Delta computation and visualization
  - Method comparison table
  - CSV export

### CSS
- `css/style.css` - Added delta visualization styles

### Documentation
- `docs/PREDICT_PHASE_B.md` (this file)

---

## Usage Guide

1. **Select T-cell type** (CD4, CD8, or Treg)
2. **Apply filters** (donor type, tissue) - click "Apply"
3. **Search for genes** using the search box
4. **Click genes** to select for knockout (max 5)
5. **Click "Apply Knockout"** to run prediction
6. **View results:**
   - Delta card shows top changes by magnitude
   - Comparison table shows all methods side-by-side
   - Full heatmaps show all 250 genes
7. **Export to CSV** for downstream analysis
