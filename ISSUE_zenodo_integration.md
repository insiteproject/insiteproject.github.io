# Issue: Integrate Zenodo Experimental Perturbation Datasets

## Problem

The current INSITE website only displays data from **scBaseCount**, which contains primarily observational disease studies (infection vs healthy controls). The perturbation filter dropdown shows limited categories:

| Current Category | Count |
|------------------|-------|
| infection | 281 |
| control | 210 |
| autoimmune | 2 |
| other | 2 |

This limits the platform's utility for researchers interested in **experimental perturbations** like:
- CRISPR knockouts
- CAR-T engineering
- Drug treatments
- Cytokine stimulation

## Available Data (Not Yet Integrated)

A Zenodo search identified **86 T-cell scRNA-seq datasets** that were evaluated and approved for inclusion:

| Category | Datasets | Examples |
|----------|----------|----------|
| Drug/Treatment | 20 | Therapeutic interventions, inhibitor studies |
| CAR-T/Engineering | 9 | CD19 CAR-T, TCR-T, engineered T cells |
| Infection | 9 | SARS-CoV-2, HIV longitudinal studies |
| Cancer/Tumor | 8 | TIL profiling, melanoma, myeloma |
| Cytokine/Stimulation | 7 | Activation assays, IFN treatment |
| Vaccination | 5 | mRNA vaccines, influenza response |
| CRISPR/Knockout | 4 | Perturbation screens, gene editing |
| Uncategorized | 23 | Various T-cell studies |

### Key Datasets for Perturbation Analysis

1. **CAR-T Studies**
   - "Single-cell expression and TCR data from CD19-specific CAR T cells in a phase I/II clinical trial"
   - Pre/post infusion comparisons available

2. **CRISPR Perturbation**
   - "Causal identification of single-cell experimental perturbation effects"
   - Direct perturbation vs control comparisons

3. **Vaccination Studies**
   - "mRNA vaccination boosts spike-specific T cell memory"
   - Pre/post vaccination timepoints

## Location of Evaluation Data

- Zenodo search results: `/fast/datawork/insite/scripts/zenodo/zenodo_search_results.json`
- Evaluation log: `/fast/datawork/insite/scripts/zenodo/zenodo_evaluation_log.md`
- Download script: `/fast/datawork/insite/scripts/zenodo/1_download_zenodo_record.py`

## Proposed Solution

### Phase 1: Download and Process Priority Datasets
1. Download top 10-15 datasets with clear perturbation categories (CAR-T, CRISPR, drug treatment)
2. Process through existing pipeline:
   - Cell typing (SingleR + marker genes)
   - HVG selection (use existing 2000-gene list)
   - UMAP computation
   - Expression profile extraction

### Phase 2: Data Integration
1. Merge processed Zenodo data with existing scBaseCount data
2. Update `datasets.json` with new samples
3. Ensure consistent gene space across all datasets

### Phase 3: Metadata Harmonization
1. Map Zenodo study metadata to perturbation categories
2. Extract timepoint information for longitudinal studies
3. Add perturbation details (target gene, drug name, dose, etc.)

## Expected Outcome

After integration, the perturbation filter should show:

| Category | Estimated Count |
|----------|-----------------|
| control | ~300 |
| infection | ~300 |
| drug_treatment | ~50-100 |
| car_t_engineering | ~30-50 |
| cytokine_stimulation | ~30-50 |
| vaccination | ~20-30 |
| genetic_perturbation | ~20-30 |
| cancer_til | ~30-50 |

## Technical Considerations

1. **Gene Space Unification**: Zenodo datasets may have different gene annotations. Need to map to common gene symbols.

2. **Batch Effects**: Consider batch correction (Harmony, scVI) when combining datasets from different sources.

3. **File Size**: Adding 50+ new datasets will increase `datasets.json` size. May need to implement lazy loading or pagination.

4. **Processing Time**: Full pipeline processing for 86 datasets could take significant compute time. Prioritize high-value datasets first.

## Priority

**High** - This significantly enhances the platform's value for T-cell engineering research, which is the core mission of INSITE.

---

*Related files:*
- Current data: `/fast/datawork/insite/website/data/datasets.json`
- Processing scripts: `/fast/datawork/insite/scripts/unification/`
- Zenodo scripts: `/fast/datawork/insite/scripts/zenodo/`
