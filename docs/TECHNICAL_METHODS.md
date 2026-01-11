# INSITE Technical Methods

## Supplementary Methods for Data Processing and Web Platform Implementation

---

## 1. Data Source

### 1.1 scBaseCount Repository

INSITE leverages the scBaseCount database, an AI agent-curated, uniformly processed repository of single-cell RNA sequencing (scRNA-seq) data. scBaseCount aggregates publicly available data from NCBI Sequence Read Archive (SRA) and Gene Expression Omnibus (GEO), providing standardized metadata extraction and consistent processing pipelines.

**Reference:** Asiaee A, et al. scBaseCount: AI agent-curated single cell data repository. bioRxiv 2025.02.27.640494 (2025).

### 1.2 T-cell Corpus Selection

The T-cell corpus was extracted from scBaseCount using the following criteria:
- Cell type annotation: T-cell or T-cell subtype (CD4+, CD8+, regulatory T-cell)
- Species: Homo sapiens
- Data format: 10x Genomics or compatible single-cell protocols
- Quality metrics: Minimum 500 genes detected per cell, maximum 20% mitochondrial content

---

## 2. Sample Selection and Quality Filtering

### 2.1 Initial Dataset

The initial T-cell corpus contained **525 samples** from various experimental studies.

### 2.2 Filtering Criteria

Samples were filtered based on completeness of metadata annotations:

1. **T-cell Subtype Filter**: Removed samples with undetermined T-cell subtype classification (labeled as "unknown")
2. **Donor Health Status Filter**: Removed samples where donor health status could not be determined from available metadata

### 2.3 Final Dataset Composition

After filtering, the final dataset contains **495 samples** with the following distribution:

| Category | Subcategory | Count |
|----------|-------------|-------|
| **T-cell Subtype** | CD4+ T-cells | 179 |
| | CD8+ T-cells | 169 |
| | Regulatory T-cells (Treg) | 147 |
| **Donor Health Status** | Healthy | 214 |
| | Diseased | 281 |

---

## 3. Donor Health Status Classification

### 3.1 Metadata Source

Donor health status was derived from the scBaseCount `disease` metadata field, which is extracted by the SRAgent AI system from SRA study descriptions, abstracts, and associated publications.

### 3.2 Classification Algorithm

A keyword-based classification approach was implemented to standardize disease annotations:

**Healthy Classification Keywords:**
- `normal`, `healthy`, `control`, `uninfected`, `mock`, `naive`

**Diseased Classification Keywords:**
- `covid`, `covid-19`, `sars-cov-2`
- `infection`, `infected`
- `disease`, `disorder`
- `cancer`, `tumor`, `carcinoma`, `lymphoma`, `leukemia`
- `hiv`, `hiv-1`
- `influenza`, `flu`
- `tuberculosis`, `tb`
- `pneumonia`
- `sepsis`
- `autoimmune`

### 3.3 Classification Logic

```
function classify_donor(disease_string):
    disease_lower = lowercase(disease_string)

    # Check healthy keywords first
    for keyword in healthy_keywords:
        if keyword in disease_lower:
            # But override if disease keyword also present
            # (e.g., "healthy control vs COVID-19 patient")
            for disease_keyword in diseased_keywords:
                if disease_keyword in disease_lower:
                    return "diseased"
            return "healthy"

    # Check disease keywords
    for keyword in diseased_keywords:
        if keyword in disease_lower:
            return "diseased"

    return "unknown"  # Excluded from final dataset
```

### 3.4 Sample-to-Metadata Matching

Metadata was linked to processed samples via SRX accession identifiers extracted from the h5ad file paths. Each processed sample path contains the SRX accession (e.g., `SRX12345678_processed.h5ad`), which was matched to the corresponding entry in the scBaseCount metadata CSV.

---

## 4. T-cell Subtype Annotation

### 4.1 Annotation Method

T-cell subtypes were determined using a two-stage approach:

1. **Reference-based Annotation (SingleR):** Cells were annotated using the SingleR algorithm with the Monaco Immune Cell Reference dataset.

2. **Marker Gene Scoring:** Classification was refined using canonical T-cell markers:

| Subtype | Positive Markers | Negative Markers |
|---------|-----------------|------------------|
| CD4+ T-cells | CD4, CD3D, CD3E | CD8A, CD8B |
| CD8+ T-cells | CD8A, CD8B, CD3D, CD3E | CD4 |
| Regulatory T-cells | FOXP3, IL2RA (CD25), CD4 | CD8A, CD8B |

### 4.2 Confidence Threshold

Cells were assigned to a subtype if:
- SingleR pruned score > 0.5, AND
- Marker gene score difference > 0.3 between top two candidates

Samples with >80% cells classified to a single subtype were assigned that subtype label.

---

## 5. Perturbation Type Categorization

### 5.1 Perturbation Categories

In the current Explore-phase release, the website **clusters unstructured scBaseCount `perturbation` text** into coarse categories for fast browsing (because the raw strings are highly variable and partially incomplete).

**Current categories (as shipped in `website/data/datasets.json`):**
- `control`
- `vaccination`
- `infection`
- `drug`
- `cytokine`
- `genetic`
- `stimulation`
- `diet_metabolic`
- `cancer`
- `autoimmune`
- `environment`
- `other`

These labels are intended to support contextual exploration (e.g., infection vs control) and will be refined as more structured perturbation datasets are integrated (see `ISSUE_zenodo_integration.md`).

**Planned finer-grained mapping (next phase):** classify experimental perturbations into standardized categories such as:

| Category | Description | Examples |
|----------|-------------|----------|
| Genetic Knockout | CRISPR/Cas9-mediated gene disruption | PDCD1 KO, LAG3 KO |
| Genetic Knockdown | RNAi/shRNA-mediated gene silencing | siTOX, shPD1 |
| Overexpression | Ectopic gene expression | CAR-T, TCR engineering |
| Small Molecule | Pharmacological treatment | Inhibitors, agonists |
| Cytokine Treatment | Cytokine stimulation/blocking | IL-2, IL-15, anti-IL-6 |
| Activation | T-cell activation protocols | anti-CD3/CD28, PMA/ionomycin |
| Infection | Pathogen exposure | Viral infection, bacterial challenge |
| Co-culture | Cell-cell interaction | Tumor co-culture, APC co-culture |
| Control | Untreated/vehicle control | Mock, DMSO, untransduced |

### 5.2 Annotation Source

Perturbation annotations were extracted from:
1. scBaseCount `perturbation` metadata field
2. Study abstracts and methods sections (via SRAgent extraction)
3. Manual curation for ambiguous cases

---

## 6. Highly Variable Gene (HVG) Selection

### 6.1 Normalization

Expression data was normalized using the standard scRNA-seq workflow:

1. **Library Size Normalization:** Counts per cell were scaled to 10,000 total counts
2. **Log Transformation:** log1p(normalized_counts)

### 6.2 HVG Selection Method

Highly variable genes were identified using the variance-stabilizing transformation (VST) method:

1. Fit a loess regression of log(variance) on log(mean) for all genes
2. Calculate standardized variance as the residual from the fitted curve
3. Rank genes by standardized variance
4. Select top **2,000 genes** as highly variable

### 6.3 Gene Filtering

Genes were pre-filtered before HVG selection:
- Minimum 3 cells expressing the gene (>0 counts)
- Excluded: mitochondrial genes (MT-*), ribosomal genes (RPS*, RPL*)
- Excluded: immunoglobulin genes (IG[HKL][VDJ]*)

---

## 7. Expression Profile Computation

### 7.1 Per-Sample Expression Profiles

For each sample, mean expression profiles were computed across the 2,000 HVGs:

```
mean_expression[gene] = mean(log1p_normalized_counts[gene, all_cells])
```

### 7.2 Aggregated Expression Profiles

When multiple samples are selected (e.g., via filters), a weighted mean expression profile is computed:

```
weighted_expression[gene] = Σ(sample_weight[i] × mean_expression[i, gene]) / Σ(sample_weight[i])

where sample_weight[i] = n_cells[i]  # Number of cells in sample i
```

### 7.3 Expression Visualization Thresholds

For heatmap visualization, expression values are binned into 5 categories based on log-normalized scRNA-seq expression ranges:

| Expression Level | Threshold Range | Visual Representation |
|-----------------|-----------------|----------------------|
| Very Low | < 0.25 | Lightest color |
| Low | 0.25 - 0.5 | Light color |
| Medium | 0.5 - 1.0 | Medium color |
| High | 1.0 - 2.0 | Dark color |
| Very High | > 2.0 | Darkest color |

These thresholds were calibrated based on the empirical distribution of expression values in the T-cell corpus, where 78.5% of values fall below 0.5 (log-normalized).

---

## 8. Dimensionality Reduction and Visualization

### 8.1 Principal Component Analysis (PCA)

PCA was performed on the HVG expression matrix:
- Centered and scaled per gene
- Top 50 principal components retained
- Explaining approximately 85% of variance

### 8.2 UMAP Embedding

Uniform Manifold Approximation and Projection (UMAP) was computed for visualization:

**Parameters:**
- n_neighbors: 30
- min_dist: 0.3
- metric: cosine
- n_components: 2

### 8.3 Per-Sample UMAP Coordinates

For web visualization, representative UMAP coordinates are stored per sample:
- Centroid of all cells in the sample
- Used for interactive scatter plot display

---

## 9. Web Platform Implementation

### 9.1 Data Format

Processed data is stored in JSON format for web delivery:

**datasets.json Structure:**
```json
{
  "id": "unique_sample_id",
  "name": "Sample Display Name",
  "srx_accession": "SRX12345678",
  "repo_id": "scbasecount",
  "n_cells": 5000,
  "tcell_type": "CD4",
  "donor_type": "healthy",
  "perturbation_type": "Control",
  "tissue": "PBMC",
  "umap_x": 1.234,
  "umap_y": -0.567,
  "mean_expression_profile": {
    "GENE1": 0.456,
    "GENE2": 1.234,
    ...
  }
}
```

### 9.2 Interactive Filtering

The web interface supports multi-dimensional filtering:
- T-cell subtype (CD4, CD8, Treg)
- Donor health status (Healthy, Diseased)
- Perturbation type (9 categories)
- Tissue of origin

Filters are applied client-side for responsive interaction.

### 9.3 Expression Heatmap

The gene expression heatmap displays:
- Top 500 genes from the 2,000 HVG list (limited for browser performance)
- Weighted mean expression across filtered samples
- Standard deviation for uncertainty visualization
- Color-coded by expression level (5-tier scale)

### 9.4 Performance Optimizations

To ensure responsive browser performance, the following limits are applied:

| Parameter | Limit | Description |
|-----------|-------|-------------|
| Display Genes | 500 | Maximum genes shown in expression heatmaps |
| UMAP Points per Sample | 100 | Maximum cells subsampled per dataset for UMAP |
| Total UMAP Points | 50,000 | Maximum total points in UMAP visualization |

These limits are configurable in the JavaScript source code (`explore.js`).

---

## 10. Software and Dependencies

### 10.1 Data Processing Pipeline

| Tool | Version | Purpose |
|------|---------|---------|
| Python | 3.10+ | Pipeline orchestration |
| Scanpy | 1.9.x | scRNA-seq analysis |
| AnnData | 0.9.x | Data structures |
| SingleR (R) | 2.0.x | Cell type annotation |
| pandas | 2.0.x | Data manipulation |

### 10.2 Web Platform

| Technology | Purpose |
|------------|---------|
| HTML5/CSS3 | Structure and styling |
| JavaScript (ES6+) | Interactive functionality |
| Bootstrap 5.3 | Responsive layout |
| Plotly.js | UMAP visualization |
| Font Awesome 6 | Icons |

---

## 11. Data Availability

### 11.1 Processed Data

The processed T-cell dataset is available through the INSITE web platform at:
- https://insiteproject.github.io

### 11.2 Source Data

Raw scRNA-seq data is available from:
- NCBI SRA/GEO (accession numbers in dataset metadata)
- scBaseCount repository: https://scbasecount.org

### 11.3 Code Availability

Processing scripts and web platform source code:
- https://github.com/insiteproject/insiteproject.github.io

---

## 12. Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | January 2025 | Initial release with 495 T-cell samples |

---

*Document prepared for INSITE: Integrating scRNA-seq for T-cell Engineering*
*DataWorks! Prize 2024 Submission*
