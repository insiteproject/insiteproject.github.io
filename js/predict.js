// Global variables
let genes = [];
let selectedGenes = new Set();
let originalExpression = {};
let originalUncertainty = {};
let predictedExpression = {};
let predictedUncertainty = {};
let datasets = [];
let currentFilters = {};
let causalMatrix = null;
let currentCellType = 'CD4';

// Load genes and initialize page
document.addEventListener('DOMContentLoaded', function() {
    // Update the existing spinner message instead of creating a new one
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
    options.forEach(option => {
        const optionElement = document.createElement('option');
        optionElement.value = option;
        optionElement.textContent = option;
        select.appendChild(optionElement);
    });
    
    // Set default selection for T-cell subtype
    if (elementId === 'predictTcellSubtype' && options.includes('CD4')) {
        select.value = 'CD4';
    }
}

// Update UMAP visualization based on filters and selected genes
async function updatePredictVisualization() {
    showSpinner('Applying filters and updating visualization...');
    
    try {
        const filters = {
            t_cell_subtype: document.getElementById('predictTcellSubtype').value,
            donor_type: document.getElementById('predictDonorType').value,
            time_point: document.getElementById('predictTimePoint').value,
            location: document.getElementById('predictTissueOrgan').value
        };

        // Store current filters for expression computation
        currentFilters = filters;
        
        // Load causal model if cell type changed
        const newCellType = filters.t_cell_subtype;
        if (newCellType !== currentCellType) {
            currentCellType = newCellType;
            await loadCausalModel(newCellType);
        }


        // Recompute expression profile when filters change
        await computeOriginalExpression();
        
        // Update selection summary
        await updatePredictSelectionSummary();
        
        // Re-render gene grid with updated expression data only if genes are loaded
        if (genes.length > 0) {
            renderGeneGrid();
        }

        // Remove UMAP image loading - will be drawn in browser instead
        // const imageName = generatePredictImageName(filters, selectedGenes, colorBy);
        // const timestamp = Date.now();
        // const imagePath = `images/predict/umaps/${imageName}?t=${timestamp}`;
        // document.getElementById('predictUmapPlot').src = imagePath;
        
    } catch (error) {
        console.error('Error updating predict visualization:', error);
    } finally {
        // Hide spinner after a minimum delay
        setTimeout(() => {
            hideSpinner();
        }, 1000);
    }
}

// Load causal model for the selected cell type
async function loadCausalModel(cellType) {
    if (!cellType || cellType === 'all') {
        causalMatrix = null;
        return;
    }
    
    try {
        // Map cell type to causal model filename
        const modelFilename = getCausalModelFilename(cellType);
        const response = await fetch(`data/${modelFilename}`);
        
        if (!response.ok) {
            console.warn(`Causal model not found for ${cellType}, using fallback`);
            causalMatrix = null;
            return;
        }
        
        const modelData = await response.json();
        causalMatrix = modelData.causal_matrix;
        console.log(`Loaded causal model for ${cellType}:`, causalMatrix.length + 'x' + causalMatrix[0].length);
        
    } catch (error) {
        console.error('Error loading causal model:', error);
        causalMatrix = null;
    }
}

// Map cell type to causal model filename
function getCausalModelFilename(cellType) {
    const cleanType = cellType.toLowerCase();
    return `${cleanType}_causal_model.json`;
}

// Generate predict image filename based on filters and selected genes
function generatePredictImageName(filters, selectedGenes, colorBy = 'tcell_subtype') {
    const parts = [];
    
    if (filters.t_cell_subtype && filters.t_cell_subtype !== 'all') {
        parts.push(`tcell-${filters.t_cell_subtype.replace(/\s+/g, '-').toLowerCase()}`);
    }
    if (filters.donor_type && filters.donor_type !== 'all') {
        parts.push(`donor-${filters.donor_type.replace(/\s+/g, '-').toLowerCase()}`);
    }
    if (filters.time_point && filters.time_point !== 'all') {
        parts.push(`time-${filters.time_point.replace(/\s+/g, '-').toLowerCase()}`);
    }
    if (filters.location && filters.location !== 'all') {
        parts.push(`tissue-${filters.location.replace(/\s+/g, '-').toLowerCase()}`);
    }
    
    // Add color by parameter
    parts.push(`color-${colorBy.replace(/_/g, '-')}`);
    
    // Add selected genes count to filename
    if (selectedGenes.size > 0) {
        parts.push(`genes-${selectedGenes.size}`);
        // Add hash of selected genes for uniqueness
        const geneHash = Array.from(selectedGenes).sort().join(',').split('').reduce((a,b) => {
            a = ((a << 5) - a) + b.charCodeAt(0);
            return a & a;
        }, 0);
        parts.push(`hash-${Math.abs(geneHash)}`);
    }
    
    return parts.length > 0 ? `${parts.join('_')}.png` : 'default.png';
}

// Load genes using DataService
async function loadGenes() {
    try {
        genes = await dataService.loadGenes();
        genes.sort(); // Sort alphabetically
        // Don't compute expression here - wait for UI to be initialized
    } catch (error) {
        console.error('Error loading genes:', error);
    }
}

// Compute expression values from filtered datasets
async function computeOriginalExpression() {
    try {
        const expressionProfile = await dataService.computeExpressionProfile(currentFilters);
        
        genes.forEach(gene => {
            if (expressionProfile[gene]) {
                originalExpression[gene] = expressionProfile[gene].mean;
                originalUncertainty[gene] = expressionProfile[gene].std;
            } else {
                // Fallback to zero if gene not found
                originalExpression[gene] = 0;
                originalUncertainty[gene] = 0;
            }
        });
        
    } catch (error) {
        console.error('Error computing original expression:', error);
        // Fallback to zero values
        genes.forEach(gene => {
            originalExpression[gene] = 0;
            originalUncertainty[gene] = 0;
        });
    }
}

// Render gene selection grid using unified heatmap function
function renderGeneGrid() {
    const container = document.getElementById('geneGrid');
    renderHeatmap(container, genes, originalExpression, originalUncertainty, {
        interactive: true,
        selectedGenes: selectedGenes,
        onGeneClick: 'toggleGene'
    });
}

// Toggle gene selection
async function toggleGene(gene) {
    if (selectedGenes.has(gene)) {
        selectedGenes.delete(gene);
    } else {
        selectedGenes.add(gene);
    }
    
    updateSelectedCount();
    renderGeneGrid();
    // Don't call updatePredictVisualization here - only update on Apply button
}

// Clear all selections
async function clearSelection() {
    selectedGenes.clear();
    updateSelectedCount();
    renderGeneGrid();
    document.getElementById('predictedCard').style.display = 'none';
    // Don't call updatePredictVisualization here - only update on Apply button
}

// Update selected gene count
function updateSelectedCount() {
    document.getElementById('selectedCount').textContent = selectedGenes.size;
}

// Update selection summary using DataService
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


// Internal version of updatePredictVisualization without spinner
async function updatePredictVisualizationInternal() {
    try {
        const filters = {
            t_cell_subtype: document.getElementById('predictTcellSubtype').value,
            donor_type: document.getElementById('predictDonorType').value,
            time_point: document.getElementById('predictTimePoint').value,
            location: document.getElementById('predictTissueOrgan').value
        };

        // Store current filters for expression computation
        currentFilters = filters;
        
        // Load causal model if cell type changed
        const newCellType = filters.t_cell_subtype;
        if (newCellType !== currentCellType) {
            currentCellType = newCellType;
            await loadCausalModel(newCellType);
        }


        // Recompute expression profile when filters change
        await computeOriginalExpression();
        
        // Update selection summary
        await updatePredictSelectionSummary();
        
        // Re-render gene grid with updated expression data only if genes are loaded
        if (genes.length > 0) {
            renderGeneGrid();
        }

        // Remove UMAP image loading - will be drawn in browser instead
        // const imageName = generatePredictImageName(filters, selectedGenes, colorBy);
        // const timestamp = Date.now();
        // const imagePath = `images/predict/umaps/${imageName}?t=${timestamp}`;
        // document.getElementById('predictUmapPlot').src = imagePath;
        
    } catch (error) {
        console.error('Error updating predict visualization:', error);
    }
}

// Apply knockout and show prediction
async function applyKnockout() {
    if (selectedGenes.size === 0) {
        alert('Please select at least one gene to knockout.');
        return;
    }
    
    showSpinner('Applying knockout and generating predictions...');
    
    try {
        generatePredictedExpression();
        renderPredictedHeatmap();
        
        // Call updatePredictVisualization without showing another spinner
        await updatePredictVisualizationInternal();
        
        document.getElementById('predictedCard').style.display = 'block';
        
    } catch (error) {
        console.error('Error applying knockout:', error);
        alert('Error applying knockout: ' + error.message);
    } finally {
        // Hide spinner after a minimum delay, then scroll
        setTimeout(() => {
            hideSpinner();
            
            // Scroll to predicted results after spinner is hidden
            setTimeout(() => {
                document.getElementById('predictedCard').scrollIntoView({ 
                    behavior: 'smooth', 
                    block: 'start' 
                });
            }, 100); // Small delay to ensure spinner fade-out completes
        }, 1000);
    }
}

// Generate predicted expression after knockout using causal model
function generatePredictedExpression() {
    if (causalMatrix && causalMatrix.length === genes.length) {
        // Use causal model for prediction
        predictedExpression = computeCausalPrediction(causalMatrix, selectedGenes, genes);
        
        // Set uncertainty to zero for knocked out genes, keep original for others
        genes.forEach(gene => {
            if (selectedGenes.has(gene)) {
                predictedUncertainty[gene] = 0;
            } else {
                predictedUncertainty[gene] = originalUncertainty[gene] || 0;
            }
        });
    } else {
        // Fallback to simple simulation if no causal model available
        genes.forEach(gene => {
            if (selectedGenes.has(gene)) {
                // Knocked out genes have zero expression and uncertainty
                predictedExpression[gene] = 0;
                predictedUncertainty[gene] = 0;
            } else {
                // For non-knocked-out genes, use original expression with some random variation
                // This simulates the downstream effects of the knockout
                const originalMean = originalExpression[gene];
                const originalStd = originalUncertainty[gene];
                
                // Add some random perturbation to simulate knockout effects
                const perturbationFactor = 0.8 + Math.random() * 0.4; // 0.8 to 1.2
                predictedExpression[gene] = Math.max(0, originalMean * perturbationFactor);
                predictedUncertainty[gene] = Math.max(0, originalStd * perturbationFactor);
            }
        });
    }
}

// Compute causal prediction using inv(I-A)*N
function computeCausalPrediction(causalMatrix, knockedOutGenes, geneList) {
    const n = geneList.length;
    
    // Create modified causal matrix A by removing outgoing edges from knocked out genes
    const A = causalMatrix.map((row, i) => {
        const gene = geneList[i];
        if (knockedOutGenes.has(gene)) {
            // Zero out this gene's outgoing edges (its row)
            return new Array(n).fill(0);
        } else {
            return [...row]; // Copy the row
        }
    });
    
    // Create identity matrix I
    const I = Array(n).fill(null).map((_, i) => 
        Array(n).fill(null).map((_, j) => i === j ? 1 : 0)
    );
    
    // Compute I - A
    const IminusA = I.map((row, i) => 
        row.map((val, j) => val - A[i][j])
    );
    
    // Generate random noise vector N
    const N = Array(n).fill(null).map(() => Math.random() * 0.5 - 0.25); // Small random values
    
    try {
        // Compute inv(I-A) * N
        const invIminusA = matrixInverse(IminusA);
        const result = matrixVectorMultiply(invIminusA, N);
        
        // Convert result to expression object and ensure non-negative values
        const prediction = {};
        geneList.forEach((gene, i) => {
            if (knockedOutGenes.has(gene)) {
                prediction[gene] = 0; // Knocked out genes have zero expression
            } else {
                // Add baseline expression and ensure non-negative
                const baselineExpression = originalExpression[gene] || 0;
                prediction[gene] = Math.max(0, baselineExpression + result[i]);
            }
        });
        
        return prediction;
        
    } catch (error) {
        console.error('Error computing causal prediction:', error);
        // Fallback to simple prediction
        const prediction = {};
        geneList.forEach(gene => {
            if (knockedOutGenes.has(gene)) {
                prediction[gene] = 0;
            } else {
                prediction[gene] = originalExpression[gene] || 0;
            }
        });
        return prediction;
    }
}

// Matrix inversion using math.js
function matrixInverse(matrix) {
    try {
        return math.inv(matrix);
    } catch (error) {
        throw new Error('Matrix is singular');
    }
}

// Matrix-vector multiplication using math.js
function matrixVectorMultiply(matrix, vector) {
    return math.multiply(matrix, vector);
}

// Render predicted expression heatmap using unified heatmap function
function renderPredictedHeatmap() {
    const container = document.getElementById('predictedHeatmap');
    renderHeatmap(container, genes, predictedExpression, predictedUncertainty, {
        interactive: false,
        selectedGenes: selectedGenes
    });
}


// Toggle heatmap view between grid and collapsed (single row)
function toggleHeatmapView(containerId) {
    let rerenderCallback;
    
    if (containerId === 'predictedHeatmap') {
        rerenderCallback = renderPredictedHeatmap;
    } else if (containerId === 'geneGrid') {
        rerenderCallback = renderGeneGrid;
    }
    
    // Use the common toggle function
    toggleHeatmapView(containerId, rerenderCallback);
}
