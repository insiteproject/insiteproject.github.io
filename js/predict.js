// Global variables
let genes = [];
let selectedGenes = new Set();
let originalExpression = {};
let originalUncertainty = {};
let predictedExpression = {};
let predictedUncertainty = {};
let datasets = [];
let currentFilters = {};

// Load genes and initialize page
document.addEventListener('DOMContentLoaded', function() {
    loadDatasets();
    loadGenes();
});

// Load datasets using DataService
async function loadDatasets() {
    try {
        datasets = await dataService.loadDatasets();
        populatePredictFilters();
        updatePredictVisualization();
    } catch (error) {
        console.error('Error loading datasets:', error);
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
}

// Update UMAP visualization based on filters and selected genes
async function updatePredictVisualization() {
    const filters = {
        t_cell_subtype: document.getElementById('predictTcellSubtype').value,
        donor_type: document.getElementById('predictDonorType').value,
        time_point: document.getElementById('predictTimePoint').value,
        location: document.getElementById('predictTissueOrgan').value
    };

    // Store current filters for expression computation
    currentFilters = filters;

    // Get color by option if it exists (only when UMAP card is visible)
    const colorByElement = document.getElementById('predictColorBy');
    const colorBy = colorByElement ? colorByElement.value : 'tcell_subtype';

    // Recompute expression profile when filters change
    await computeOriginalExpression();
    
    // Update selection summary
    await updatePredictSelectionSummary();
    
    // Re-render gene grid with updated expression data
    renderGeneGrid();

    // Remove UMAP image loading - will be drawn in browser instead
    // const imageName = generatePredictImageName(filters, selectedGenes, colorBy);
    // const timestamp = Date.now();
    // const imagePath = `images/predict/umaps/${imageName}?t=${timestamp}`;
    // document.getElementById('predictUmapPlot').src = imagePath;
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
        // Initialize with default filters
        currentFilters = {
            t_cell_subtype: 'all',
            donor_type: 'all', 
            time_point: 'all',
            location: 'all'
        };
        await computeOriginalExpression();
        renderGeneGrid();
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
    renderHeatmap(container, originalExpression, true);
}

// Get CSS class based on expression level
function getExpressionClass(expression) {
    if (expression < 2) return 'expr-low';
    if (expression < 4) return 'expr-low-med';
    if (expression < 6) return 'expr-medium';
    if (expression < 8) return 'expr-med-high';
    return 'expr-high';
}

// Get CSS class based on uncertainty level (using different colormap)
function getUncertaintyClass(uncertainty) {
    if (uncertainty < 0.4) return 'uncert-very-low';
    if (uncertainty < 0.8) return 'uncert-low';
    if (uncertainty < 1.2) return 'uncert-medium';
    if (uncertainty < 1.6) return 'uncert-high';
    return 'uncert-very-high';
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
    await updatePredictVisualization(); // Update UMAP when genes change
}

// Clear all selections
async function clearSelection() {
    selectedGenes.clear();
    updateSelectedCount();
    renderGeneGrid();
    document.getElementById('predictedCard').style.display = 'none';
    document.getElementById('predictUmapCard').style.display = 'none';
    await updatePredictVisualization();
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

// Handle color by change event (wrapper for async function)
async function handleColorByChange() {
    await updatePredictVisualization();
}

// Apply knockout and show prediction
async function applyKnockout() {
    if (selectedGenes.size === 0) {
        alert('Please select at least one gene to knockout.');
        return;
    }
    
    generatePredictedExpression();
    renderPredictedHeatmap();
    await updatePredictVisualization(); // Update UMAP with prediction
    document.getElementById('predictedCard').style.display = 'block';
    document.getElementById('predictUmapCard').style.display = 'block';
    
    // Scroll to predicted results
    document.getElementById('predictedCard').scrollIntoView({ 
        behavior: 'smooth', 
        block: 'start' 
    });
}

// Generate predicted expression after knockout
function generatePredictedExpression() {
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

// Render predicted expression heatmap using unified heatmap function
function renderPredictedHeatmap() {
    const container = document.getElementById('predictedHeatmap');
    renderHeatmap(container, predictedExpression, false);
}

// Generic heatmap renderer
function renderHeatmap(container, expressionData, interactive) {
    // Check if container is in collapsed state
    const isCollapsed = container.classList.contains('heatmap-collapsed');
    const gridClass = isCollapsed ? 'heatmap-grid-collapsed' : 'heatmap-grid';
    
    // Determine which uncertainty data to use
    const uncertaintyData = expressionData === originalExpression ? originalUncertainty : predictedUncertainty;
    
    let html = generateColorbarHTML() + `<div class="${gridClass}">`;
    
    genes.forEach(gene => {
        const expression = expressionData[gene];
        const uncertainty = uncertaintyData[gene];
        const expressionClass = getExpressionClass(expression);
        const uncertaintyClass = getUncertaintyClass(uncertainty);
        const isKnockedOut = selectedGenes.has(gene);
        const interactiveClass = interactive ? 'interactive' : '';
        const clickHandler = interactive ? `onclick="toggleGene('${gene}')"` : '';
        
        html += `
            <div class="heatmap-cell ${expressionClass} ${interactiveClass} ${isKnockedOut ? 'knocked-out' : ''}" 
                 ${clickHandler}
                 data-gene="${gene}"
                 title="${gene}: Mean=${expression.toFixed(2)}, Std=${uncertainty.toFixed(2)}">
                <div class="heatmap-gene-name">${gene}</div>
                <div class="heatmap-expression">${expression.toFixed(1)}</div>
                <div class="uncertainty-bar ${uncertaintyClass}">
                    <span class="uncertainty-value">${uncertainty.toFixed(1)}</span>
                </div>
                ${isKnockedOut ? '<div class="knockout-slash"></div>' : ''}
            </div>
        `;
    });
    
    html += '</div>';
    container.innerHTML = html;
}

// Toggle heatmap view between grid and collapsed (single row)
function toggleHeatmapView(containerId) {
    const container = document.getElementById(containerId);
    let button;
    
    // Handle different button ID patterns
    if (containerId === 'geneGrid') {
        button = document.getElementById('geneGridCollapseBtn');
    } else {
        button = document.getElementById(containerId.replace('Heatmap', 'CollapseBtn'));
    }
    
    const isCollapsed = container.classList.contains('heatmap-collapsed');
    
    if (isCollapsed) {
        // Expand to grid view
        container.classList.remove('heatmap-collapsed');
        button.innerHTML = '<i class="fas fa-compress-alt"></i> Collapse';
    } else {
        // Collapse to single row
        container.classList.add('heatmap-collapsed');
        button.innerHTML = '<i class="fas fa-expand-alt"></i> Expand';
    }
    
    // Re-render the appropriate visualization
    if (containerId === 'predictedHeatmap') {
        renderPredictedHeatmap();
    } else if (containerId === 'geneGrid') {
        renderGeneGrid();
    }
}
