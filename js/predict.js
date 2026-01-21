// Global variables
let genes = [];
let selectedGenes = new Set();
let originalExpression = {};
let originalUncertainty = {};
let datasets = [];
let currentFilters = {};
let currentCellType = 'CD4';

// Phase B: Model storage
const models = {
    gaussian: null,
    invariant: null,
    dae: null
};

// Phase B: Baseline caching
const baselineCache = new Map();

// Phase B: Prediction results storage
let predictions = {
    naive: {},
    gaussian: {},
    invariant: {},
    dae: {}
};

// Load genes and initialize page
document.addEventListener('DOMContentLoaded', function() {
    const existingSpinner = document.getElementById('loadingSpinner');
    if (existingSpinner) {
        const messageElement = existingSpinner.querySelector('.loading-message');
        if (messageElement) {
            messageElement.textContent = 'Loading datasets and genes...';
        }
    }
    loadDatasets();
    loadGenes();
});

// Load datasets using DataService
async function loadDatasets() {
    try {
        datasets = await dataService.loadDatasets();
        await populatePredictFilters();
        await updatePredictVisualization();
    } catch (error) {
        console.error('Error loading datasets:', error);
    } finally {
        hideSpinner();
    }
}

// Populate filter dropdowns using DataService
async function populatePredictFilters() {
    try {
        const tcellSubtypes = await dataService.getUniqueValues('t_cell_subtype');
        const donorTypes = await dataService.getUniqueValues('donor_type');
        const timePoints = await dataService.getUniqueValues('time_point');
        const tissueOrgans = await dataService.getUniqueValues('location');

        populatePredictSelect('predictTcellSubtype', tcellSubtypes);
        populatePredictSelect('predictDonorType', donorTypes);
        populatePredictSelect('predictTimePoint', timePoints);
        populatePredictSelect('predictTissueOrgan', tissueOrgans);
    } catch (error) {
        console.error('Error populating predict filters:', error);
    }
}

// Helper function to populate select elements
function populatePredictSelect(elementId, options) {
    const select = document.getElementById(elementId);
    if (!select) return;
    options.forEach(option => {
        const optionElement = document.createElement('option');
        optionElement.value = option;
        optionElement.textContent = option;
        select.appendChild(optionElement);
    });

    if (elementId === 'predictTcellSubtype' && options.includes('CD4')) {
        select.value = 'CD4';
    }
}

// Update visualization based on filters
async function updatePredictVisualization() {
    showSpinner('Applying filters and updating visualization...');

    try {
        const filters = {
            t_cell_subtype: document.getElementById('predictTcellSubtype').value,
            donor_type: document.getElementById('predictDonorType').value,
            time_point: document.getElementById('predictTimePoint').value,
            location: document.getElementById('predictTissueOrgan').value
        };

        currentFilters = filters;

        // Load models if cell type changed
        const newCellType = filters.t_cell_subtype;
        if (newCellType !== currentCellType && newCellType !== 'all') {
            currentCellType = newCellType;
            await loadAllModels(newCellType);
        }

        // Compute baseline expression
        await computeOriginalExpression();

        // Update selection summary
        await updatePredictSelectionSummary();

        // Re-render gene grid
        if (genes.length > 0) {
            renderGeneGrid();
        }

    } catch (error) {
        console.error('Error updating predict visualization:', error);
    } finally {
        setTimeout(() => hideSpinner(), 500);
    }
}

// Phase B: Load all prediction models for a cell type
async function loadAllModels(cellType) {
    if (!cellType || cellType === 'all') {
        models.gaussian = null;
        models.invariant = null;
        models.dae = null;
        return;
    }

    const ct = cellType.toLowerCase();
    const baseUrl = `${dataService.baseUrl}data/predict/models/`;

    try {
        // Load all three models in parallel
        const [gaussianRes, invariantRes, daeRes] = await Promise.all([
            fetch(`${baseUrl}gaussian_conditional_${ct}.json`),
            fetch(`${baseUrl}invariant_within_bioproject_${ct}.json`),
            fetch(`${baseUrl}dae_denoiser_${ct}.json`)
        ]);

        if (gaussianRes.ok) {
            models.gaussian = await gaussianRes.json();
            console.log(`Loaded Gaussian model for ${cellType}`);
        }
        if (invariantRes.ok) {
            models.invariant = await invariantRes.json();
            console.log(`Loaded Invariant model for ${cellType}`);
        }
        if (daeRes.ok) {
            models.dae = await daeRes.json();
            console.log(`Loaded DAE model for ${cellType}`);
        }

    } catch (error) {
        console.error('Error loading models:', error);
    }
}

// Load genes using DataService
async function loadGenes() {
    try {
        genes = await dataService.loadGenes();
    } catch (error) {
        console.error('Error loading genes:', error);
    }
}

// Phase B: Get filter key for caching
function getFilterKey(filters) {
    return JSON.stringify({
        t_cell_subtype: filters.t_cell_subtype,
        donor_type: filters.donor_type,
        time_point: filters.time_point,
        location: filters.location
    });
}

// Compute expression values from filtered datasets (with caching)
async function computeOriginalExpression() {
    const filterKey = getFilterKey(currentFilters);

    // Check cache first
    if (baselineCache.has(filterKey)) {
        const cached = baselineCache.get(filterKey);
        originalExpression = cached.expression;
        originalUncertainty = cached.uncertainty;
        return;
    }

    try {
        const expressionProfile = await dataService.computeExpressionProfile(currentFilters);

        genes.forEach(gene => {
            if (expressionProfile[gene]) {
                originalExpression[gene] = expressionProfile[gene].mean;
                originalUncertainty[gene] = expressionProfile[gene].std;
            } else {
                originalExpression[gene] = 0;
                originalUncertainty[gene] = 0;
            }
        });

        // Cache the result
        baselineCache.set(filterKey, {
            expression: { ...originalExpression },
            uncertainty: { ...originalUncertainty }
        });

    } catch (error) {
        console.error('Error computing original expression:', error);
        genes.forEach(gene => {
            originalExpression[gene] = 0;
            originalUncertainty[gene] = 0;
        });
    }
}

// Render gene selection grid
function renderGeneGrid() {
    const container = document.getElementById('geneGrid');
    renderHeatmap(container, genes, originalExpression, originalUncertainty, {
        interactive: true,
        selectedGenes: selectedGenes,
        onGeneClick: 'toggleGene'
    });
}

// Phase B: Filter gene grid by search term
function filterGeneGrid(searchTerm) {
    const container = document.getElementById('geneGrid');
    const cells = container.querySelectorAll('.heatmap-cell');
    const term = searchTerm.toLowerCase().trim();

    cells.forEach(cell => {
        const geneName = cell.getAttribute('data-gene') || cell.textContent.trim();
        if (term === '' || geneName.toLowerCase().includes(term)) {
            cell.classList.remove('gene-hidden');
        } else {
            cell.classList.add('gene-hidden');
        }
    });
}

// Toggle gene selection
async function toggleGene(gene) {
    if (selectedGenes.has(gene)) {
        selectedGenes.delete(gene);
    } else {
        if (selectedGenes.size >= 5) {
            alert('Maximum 5 genes can be selected for knockout.');
            return;
        }
        selectedGenes.add(gene);
    }

    updateSelectedCount();
    renderGeneGrid();
}

// Clear all selections
async function clearSelection() {
    selectedGenes.clear();
    updateSelectedCount();
    renderGeneGrid();
    document.getElementById('predictedCard').style.display = 'none';
    document.getElementById('deltaCard').style.display = 'none';
    document.getElementById('comparisonCard').style.display = 'none';
}

// Update selected gene count
function updateSelectedCount() {
    document.getElementById('selectedCount').textContent = selectedGenes.size;
}

// Update selection summary
async function updatePredictSelectionSummary() {
    try {
        const filteredDatasets = await dataService.filterDatasets(currentFilters);
        const stats = await dataService.getStatistics(filteredDatasets);
        document.getElementById('predictDatasetCount').textContent = stats.recordCount;
        document.getElementById('predictCellCount').textContent = stats.totalCells.toLocaleString();
    } catch (error) {
        console.error('Error updating predict selection summary:', error);
        document.getElementById('predictDatasetCount').textContent = '0';
        document.getElementById('predictCellCount').textContent = '0';
    }
}

// Phase B: Apply knockout and run all prediction methods
async function applyKnockout() {
    if (selectedGenes.size === 0) {
        alert('Please select at least one gene to knockout.');
        return;
    }

    showSpinner('Running predictions with all methods...');

    try {
        // Run all 4 prediction methods
        predictions.naive = predictNaive();
        predictions.gaussian = predictGaussian();
        predictions.invariant = predictInvariant();
        predictions.dae = predictDAE();

        // Render results
        renderPredictedHeatmap();
        renderDeltaVisualization();
        renderComparisonTable();

        // Show result cards
        document.getElementById('predictedCard').style.display = 'block';
        document.getElementById('deltaCard').style.display = 'block';
        document.getElementById('comparisonCard').style.display = 'block';

    } catch (error) {
        console.error('Error applying knockout:', error);
        alert('Error applying knockout: ' + error.message);
    } finally {
        setTimeout(() => {
            hideSpinner();
            document.getElementById('predictedCard').scrollIntoView({
                behavior: 'smooth',
                block: 'start'
            });
        }, 500);
    }
}

// Phase B: Method 1 - Naive KO (just set to zero)
function predictNaive() {
    const result = {};
    genes.forEach(gene => {
        if (selectedGenes.has(gene)) {
            result[gene] = 0;
        } else {
            result[gene] = originalExpression[gene] || 0;
        }
    });
    return result;
}

// Phase B: Method 2 - Gaussian Conditional
function predictGaussian() {
    if (!models.gaussian || !models.gaussian.precision) {
        console.warn('Gaussian model not loaded, falling back to naive');
        return predictNaive();
    }
    return predictWithPrecision(models.gaussian.precision);
}

// Phase B: Method 3 - Invariant Precision
function predictInvariant() {
    if (!models.invariant || !models.invariant.precision) {
        console.warn('Invariant model not loaded, falling back to naive');
        return predictNaive();
    }
    return predictWithPrecision(models.invariant.precision);
}

// Phase B: Precision-based prediction (Gaussian conditional)
function predictWithPrecision(precision) {
    const n = genes.length;
    const koIndices = [];
    const remainingIndices = [];

    genes.forEach((gene, i) => {
        if (selectedGenes.has(gene)) {
            koIndices.push(i);
        } else {
            remainingIndices.push(i);
        }
    });

    if (koIndices.length === 0) return { ...originalExpression };

    // Convert to log space
    const mu = genes.map(g => Math.log1p(originalExpression[g] || 0));

    // Extract P_SS (KO × KO submatrix)
    const k = koIndices.length;
    const Pss = [];
    for (let i = 0; i < k; i++) {
        Pss[i] = [];
        for (let j = 0; j < k; j++) {
            Pss[i][j] = precision[koIndices[i]][koIndices[j]];
        }
    }

    // Add ridge for stability
    const ridge = 1e-6;
    for (let i = 0; i < k; i++) {
        Pss[i][i] += ridge;
    }

    // RHS: (0 - mu_S)
    const rhs = koIndices.map(i => 0 - mu[i]);

    // Solve Pss * alpha = rhs
    const alpha = solveLinearSystem(Pss, rhs);

    // Compute predictions for remaining genes
    const predLog = [...mu];

    // For KO genes: set to 0 in log space
    koIndices.forEach(i => {
        predLog[i] = 0;
    });

    // For remaining genes: mu_R - P_RS @ alpha
    remainingIndices.forEach(ri => {
        let adjustment = 0;
        for (let j = 0; j < k; j++) {
            adjustment += precision[ri][koIndices[j]] * alpha[j];
        }
        predLog[ri] = mu[ri] - adjustment;
    });

    // Convert back to counts
    const result = {};
    genes.forEach((gene, i) => {
        result[gene] = Math.max(0, Math.expm1(predLog[i]));
    });

    return result;
}

// Phase B: Solve linear system using Gaussian elimination
function solveLinearSystem(A, b) {
    const n = A.length;
    const M = A.map((row, i) => [...row, b[i]]);

    // Forward elimination
    for (let i = 0; i < n; i++) {
        // Find pivot
        let maxRow = i;
        for (let k = i + 1; k < n; k++) {
            if (Math.abs(M[k][i]) > Math.abs(M[maxRow][i])) {
                maxRow = k;
            }
        }
        [M[i], M[maxRow]] = [M[maxRow], M[i]];

        // Eliminate
        for (let k = i + 1; k < n; k++) {
            const factor = M[k][i] / M[i][i];
            for (let j = i; j <= n; j++) {
                M[k][j] -= factor * M[i][j];
            }
        }
    }

    // Back substitution
    const x = new Array(n).fill(0);
    for (let i = n - 1; i >= 0; i--) {
        x[i] = M[i][n];
        for (let j = i + 1; j < n; j++) {
            x[i] -= M[i][j] * x[j];
        }
        x[i] /= M[i][i];
    }

    return x;
}

// Phase B: Method 4 - DAE (Denoising Autoencoder)
function predictDAE() {
    if (!models.dae || !models.dae.mlp) {
        console.warn('DAE model not loaded, falling back to naive');
        return predictNaive();
    }

    // Convert to log space and create input
    let x = genes.map(g => Math.log1p(originalExpression[g] || 0));

    // Pre-clamp KO genes to 0
    genes.forEach((gene, i) => {
        if (selectedGenes.has(gene)) {
            x[i] = 0;
        }
    });

    // Forward pass through MLP
    const layers = models.dae.mlp.layers;
    for (const layer of layers) {
        x = mlpLayerForward(x, layer.weight, layer.bias, layer.activation);
    }

    // Post-clamp KO genes to 0
    genes.forEach((gene, i) => {
        if (selectedGenes.has(gene)) {
            x[i] = 0;
        }
    });

    // Convert back to counts
    const result = {};
    genes.forEach((gene, i) => {
        result[gene] = Math.max(0, Math.expm1(x[i]));
    });

    return result;
}

// Phase B: MLP layer forward pass
function mlpLayerForward(input, weight, bias, activation) {
    const outputSize = weight.length;
    const output = new Array(outputSize).fill(0);

    for (let i = 0; i < outputSize; i++) {
        output[i] = bias[i];
        for (let j = 0; j < input.length; j++) {
            output[i] += weight[i][j] * input[j];
        }

        // Apply activation
        if (activation === 'gelu') {
            output[i] = gelu(output[i]);
        }
        // 'linear' or undefined: no activation
    }

    return output;
}

// Phase B: GELU activation function
function gelu(x) {
    return 0.5 * x * (1 + Math.tanh(Math.sqrt(2 / Math.PI) * (x + 0.044715 * x * x * x)));
}

// Phase B: Render predicted heatmap (using Gaussian as default display)
function renderPredictedHeatmap() {
    const container = document.getElementById('predictedHeatmap');
    const method = document.getElementById('displayMethod')?.value || 'gaussian';
    const expression = predictions[method] || predictions.gaussian || predictions.naive;

    // Compute uncertainty (0 for KO genes)
    const uncertainty = {};
    genes.forEach(gene => {
        uncertainty[gene] = selectedGenes.has(gene) ? 0 : (originalUncertainty[gene] || 0);
    });

    renderHeatmap(container, genes, expression, uncertainty, {
        interactive: false,
        selectedGenes: selectedGenes
    });
}

// Phase B: Render delta visualization
function renderDeltaVisualization() {
    const method = document.getElementById('displayMethod')?.value || 'gaussian';
    const pred = predictions[method] || predictions.gaussian;

    // Calculate deltas
    const deltas = [];
    genes.forEach(gene => {
        const delta = (pred[gene] || 0) - (originalExpression[gene] || 0);
        deltas.push({ gene, delta, isKO: selectedGenes.has(gene) });
    });

    // Sort by delta (excluding KO genes for ranking)
    const nonKO = deltas.filter(d => !d.isKO);
    nonKO.sort((a, b) => b.delta - a.delta);

    // Top upregulated
    const topUp = nonKO.slice(0, 10);
    renderDeltaTable('topUpregulated', topUp, 'up');

    // Top downregulated
    const topDown = nonKO.slice(-10).reverse();
    renderDeltaTable('topDownregulated', topDown, 'down');
}

// Phase B: Render delta table
function renderDeltaTable(containerId, items, direction) {
    const container = document.getElementById(containerId);
    if (!container) return;

    let html = '<table class="table table-sm table-bordered mb-0">';
    html += '<thead class="table-light"><tr><th>Gene</th><th>Δ Expression</th></tr></thead>';
    html += '<tbody>';

    items.forEach(item => {
        const colorClass = item.delta > 0 ? 'text-success' : 'text-danger';
        const sign = item.delta > 0 ? '+' : '';
        html += `<tr>
            <td><strong>${item.gene}</strong></td>
            <td class="${colorClass}">${sign}${item.delta.toFixed(3)}</td>
        </tr>`;
    });

    html += '</tbody></table>';
    container.innerHTML = html;
}

// Phase B: Render method comparison table
function renderComparisonTable() {
    const container = document.getElementById('comparisonTable');
    if (!container) return;

    // Get top changed genes (by absolute delta from Gaussian)
    const gaussianDeltas = genes.map(gene => ({
        gene,
        delta: Math.abs((predictions.gaussian[gene] || 0) - (originalExpression[gene] || 0)),
        isKO: selectedGenes.has(gene)
    }));

    gaussianDeltas.sort((a, b) => b.delta - a.delta);
    const topGenes = gaussianDeltas.slice(0, 20).map(d => d.gene);

    let html = '<table class="table table-sm table-bordered">';
    html += '<thead class="table-light"><tr>';
    html += '<th>Gene</th><th>Baseline</th>';
    html += '<th>Naive</th><th>Gaussian</th><th>Invariant</th><th>DAE</th>';
    html += '</tr></thead><tbody>';

    topGenes.forEach(gene => {
        const baseline = originalExpression[gene] || 0;
        const isKO = selectedGenes.has(gene);
        const rowClass = isKO ? 'table-warning' : '';

        html += `<tr class="${rowClass}">`;
        html += `<td><strong>${gene}</strong>${isKO ? ' <span class="badge bg-danger">KO</span>' : ''}</td>`;
        html += `<td>${baseline.toFixed(2)}</td>`;

        ['naive', 'gaussian', 'invariant', 'dae'].forEach(method => {
            const pred = predictions[method][gene] || 0;
            const delta = pred - baseline;
            const colorClass = delta > 0.01 ? 'text-success' : (delta < -0.01 ? 'text-danger' : '');
            html += `<td class="${colorClass}">${pred.toFixed(2)}</td>`;
        });

        html += '</tr>';
    });

    html += '</tbody></table>';
    container.innerHTML = html;
}

// Phase B: Export predictions to CSV
function exportPredictions() {
    if (selectedGenes.size === 0) {
        alert('Please run a prediction first.');
        return;
    }

    let csv = 'Gene,Baseline,Naive,Gaussian,Invariant,DAE,Delta_Naive,Delta_Gaussian,Delta_Invariant,Delta_DAE,Is_KO\n';

    genes.forEach(gene => {
        const baseline = originalExpression[gene] || 0;
        const naive = predictions.naive[gene] || 0;
        const gaussian = predictions.gaussian[gene] || 0;
        const invariant = predictions.invariant[gene] || 0;
        const dae = predictions.dae[gene] || 0;
        const isKO = selectedGenes.has(gene) ? 1 : 0;

        csv += `${gene},${baseline},${naive},${gaussian},${invariant},${dae},`;
        csv += `${naive - baseline},${gaussian - baseline},${invariant - baseline},${dae - baseline},${isKO}\n`;
    });

    // Download
    const blob = new Blob([csv], { type: 'text/csv' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `knockout_predictions_${currentCellType}_${Array.from(selectedGenes).join('_')}.csv`;
    a.click();
    URL.revokeObjectURL(url);
}

// Toggle heatmap view
function togglePredictHeatmapView(containerId) {
    let rerenderCallback;

    if (containerId === 'predictedHeatmap') {
        rerenderCallback = renderPredictedHeatmap;
    } else if (containerId === 'geneGrid') {
        rerenderCallback = renderGeneGrid;
    }

    toggleHeatmapView(containerId, rerenderCallback);
}

// Change display method for predicted heatmap
function changeDisplayMethod() {
    renderPredictedHeatmap();
    renderDeltaVisualization();
}
