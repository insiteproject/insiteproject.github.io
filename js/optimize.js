// Global variables
let datasets = [];
let genes = [];
let selectedTargetGenes = new Set();
let knockoutBudget = 2;
let currentFilters = {};
let originalExpression = {};
let originalUncertainty = {};
let editedExpression = {}; // Track user edits to expression values

// Load data and initialize page
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
    updateBudgetDisplay();
});

// Load datasets using DataService
async function loadDatasets() {
    try {
        datasets = await dataService.loadDatasets();
        populateOptimizeFilters();
        updateOptimizeVisualization();
    } catch (error) {
        console.error('Error loading datasets:', error);
    } finally {
        hideSpinner();
    }
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

// Populate filter dropdowns using DataService
async function populateOptimizeFilters() {
    try {
        const tcellSubtypes = await dataService.getUniqueValues('t_cell_subtype');
        const donorTypes = await dataService.getUniqueValues('donor_type');
        const timePoints = await dataService.getUniqueValues('time_point');
        const tissueOrgans = await dataService.getUniqueValues('location');

        populateOptimizeSelect('optimizeTcellSubtype', tcellSubtypes);
        populateOptimizeSelect('optimizeDonorType', donorTypes);
        populateOptimizeSelect('optimizeTimePoint', timePoints);
        populateOptimizeSelect('optimizeTissueOrgan', tissueOrgans);
    } catch (error) {
        console.error('Error populating optimize filters:', error);
    }
}

// Helper function to populate select elements
function populateOptimizeSelect(elementId, options) {
    const select = document.getElementById(elementId);
    options.forEach(option => {
        const optionElement = document.createElement('option');
        optionElement.value = option;
        optionElement.textContent = option;
        select.appendChild(optionElement);
    });
    
    // Set default selection for T-cell subtype
    if (elementId === 'optimizeTcellSubtype' && options.includes('CD4')) {
        select.value = 'CD4';
    }
}

// Update selection summary using DataService
async function updateOptimizeSelectionSummary(filters) {
    try {
        const filteredDatasets = await dataService.filterDatasets(filters);
        const stats = await dataService.getStatistics(filteredDatasets);
        document.getElementById('optimizeDatasetCount').textContent = stats.recordCount;
        document.getElementById('optimizeCellCount').textContent = stats.totalCells.toLocaleString();
    } catch (error) {
        console.error('Error updating optimize selection summary:', error);
        document.getElementById('optimizeDatasetCount').textContent = '0';
        document.getElementById('optimizeCellCount').textContent = '0';
    }
}

// Update visualization based on filters
async function updateOptimizeVisualization() {
    showSpinner('Applying filters and updating visualization...');
    
    try {
        const filters = {
            t_cell_subtype: document.getElementById('optimizeTcellSubtype').value,
            donor_type: document.getElementById('optimizeDonorType').value,
            time_point: document.getElementById('optimizeTimePoint').value,
            location: document.getElementById('optimizeTissueOrgan').value
        };

        // Store current filters for expression computation
        currentFilters = filters;

        // Update selection summary
        await updateOptimizeSelectionSummary(filters);

        // Recompute expression profile when filters change
        await computeOriginalExpression();
        
        // Re-render gene heatmap with updated expression data only if genes are loaded
        if (genes.length > 0) {
            renderOptimizeHeatmap();
        }
        
        // Reset target selection when filters change
        resetTargetSelection();
        
    } catch (error) {
        console.error('Error updating optimize visualization:', error);
    } finally {
        // Hide spinner after a minimum delay
        setTimeout(() => {
            hideSpinner();
        }, 1000);
    }
}


// Update budget display
function updateBudgetDisplay() {
    const input = document.getElementById('knockoutBudget');
    let value = parseInt(input.value);
    
    // Enforce max value of 10
    if (value > 10) {
        value = 10;
        input.value = 10;
    }
    
    knockoutBudget = value;
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
            // Initialize edited expression to match original
            editedExpression[gene] = originalExpression[gene];
        });
        
    } catch (error) {
        console.error('Error computing original expression:', error);
        // Fallback to zero values
        genes.forEach(gene => {
            originalExpression[gene] = 0;
            originalUncertainty[gene] = 0;
            editedExpression[gene] = 0;
        });
    }
}

// Render gene expression heatmap with editable controls
function renderOptimizeHeatmap() {
    const container = document.getElementById('optimizeHeatmap');
    
    // Check if container is in collapsed state
    const isCollapsed = container.classList.contains('heatmap-collapsed');
    const gridClass = isCollapsed ? 'heatmap-grid-collapsed' : 'heatmap-grid';
    
    let html = generateColorbarHTML() + `<div class="${gridClass}">`;
    
    genes.forEach(gene => {
        const expression = editedExpression[gene] || 0;
        const uncertainty = originalUncertainty[gene] || 0;
        const expressionClass = getExpressionClass(expression);
        const uncertaintyClass = getUncertaintyClass(uncertainty);
        const isSelected = selectedTargetGenes.has(gene);
        const isEdited = Math.abs(editedExpression[gene] - originalExpression[gene]) > 0.01;
        
        html += `
            <div class="heatmap-cell editable-cell ${expressionClass} ${isSelected ? 'knocked-out' : ''} ${isEdited ? 'edited-cell' : ''}" 
                 data-gene="${gene}"
                 data-expression="${expression.toFixed(2)}"
                 data-uncertainty="${uncertainty.toFixed(2)}"
                 onmouseenter="showHeatmapTooltip(event)"
                 onmouseleave="hideHeatmapTooltip()">
                <div class="heatmap-gene-name">${gene}</div>
                <div class="uncertainty-dot ${uncertaintyClass}"></div>
                <div class="expression-controls">
                    <button class="expression-btn decrease-btn" onclick="adjustExpression('${gene}', -0.5)" title="Decrease expression">-</button>
                    <span class="expression-value" onclick="toggleTargetGene('${gene}')" title="Click to select/deselect">${expression.toFixed(1)}</span>
                    <button class="expression-btn increase-btn" onclick="adjustExpression('${gene}', 0.5)" title="Increase expression">+</button>
                </div>
                ${isSelected ? '<div class="knockout-slash"></div>' : ''}
                ${isEdited ? '<div class="edited-indicator">✎</div>' : ''}
            </div>
        `;
    });
    
    html += '</div>';
    container.innerHTML = html;
}

// Toggle gene selection for targeting
async function toggleTargetGene(gene) {
    if (selectedTargetGenes.has(gene)) {
        selectedTargetGenes.delete(gene);
    } else {
        selectedTargetGenes.add(gene);
    }
    
    updateSelectedTargetCount();
    renderOptimizeHeatmap();
}

// Update selected target gene count
function updateSelectedTargetCount() {
    document.getElementById('selectedTargetCount').textContent = selectedTargetGenes.size;
}

// Reset target selection
function resetTargetSelection() {
    selectedTargetGenes.clear();
    updateSelectedTargetCount();
    renderOptimizeHeatmap();
}

// Adjust gene expression value
function adjustExpression(gene, delta) {
    const currentValue = editedExpression[gene] || 0;
    const newValue = Math.max(0, Math.min(10, currentValue + delta)); // Clamp between 0 and 10
    editedExpression[gene] = newValue;
    
    // Re-render the heatmap to show updated values
    renderOptimizeHeatmap();
}

// Reset all expression values to original
function resetExpressionValues() {
    genes.forEach(gene => {
        editedExpression[gene] = originalExpression[gene];
    });
    renderOptimizeHeatmap();
}

// Run optimization algorithm
function runOptimization() {
    // Check if any genes have been selected as targets
    const hasSelectedTargets = selectedTargetGenes.size > 0;
    
    // Check if any expression values have been edited
    const hasEditedExpression = genes.some(gene => 
        Math.abs(editedExpression[gene] - originalExpression[gene]) > 0.01
    );
    
    if (!hasSelectedTargets && !hasEditedExpression) {
        alert('Please select target genes or modify expression values to run optimization.');
        return;
    }
    
    showSpinner('Running optimization algorithm...');
    
    try {
        // Simulate optimization process
        setTimeout(() => {
            const results = generateOptimizationResults();
            displayResults(results);
            
            // Hide spinner and then scroll
            setTimeout(() => {
                hideSpinner();
                
                // Scroll to results after spinner is hidden
                setTimeout(() => {
                    document.getElementById('resultsCard').scrollIntoView({ 
                        behavior: 'smooth', 
                        block: 'start' 
                    });
                }, 100); // Small delay to ensure spinner fade-out completes
            }, 1000);
        }, 2000); // 2 second delay to simulate computation
        
    } catch (error) {
        console.error('Error running optimization:', error);
        alert('Error running optimization: ' + error.message);
        hideSpinner();
    }
}

// Generate mock optimization results
function generateOptimizationResults() {
    const results = [];
    
    // Generate 5 different gene combinations
    for (let i = 0; i < 5; i++) {
        const geneCombination = [];
        const usedGenes = new Set();
        
        // Select random genes for this combination (based on budget)
        while (geneCombination.length < knockoutBudget) {
            const randomGene = genes[Math.floor(Math.random() * genes.length)];
            if (!usedGenes.has(randomGene)) {
                geneCombination.push(randomGene);
                usedGenes.add(randomGene);
            }
        }
        
        // Generate realistic scores (higher for first results)
        const baseScore = 95 - (i * 8) + Math.random() * 5;
        const score = Math.max(60, Math.min(100, baseScore));
        
        // Generate confidence (correlated with score but with some variation)
        const confidence = Math.max(70, Math.min(99, score - 5 + Math.random() * 10));
        
        results.push({
            rank: i + 1,
            genes: geneCombination,
            score: score.toFixed(1),
            confidence: confidence.toFixed(0)
        });
    }
    
    return results;
}

// Display optimization results
function displayResults(results) {
    const tableBody = document.getElementById('resultsTableBody');
    
    // Color scheme for results
    const resultColors = ['#ff6b6b', '#4ecdc4', '#45b7d1', '#96ceb4', '#feca57'];
    
    const html = results.map((result, index) => {
        const scoreClass = parseFloat(result.score) > 85 ? 'text-success' : 
                          parseFloat(result.score) > 75 ? 'text-warning' : 'text-danger';
        const confidenceClass = parseInt(result.confidence) > 90 ? 'text-success' : 
                               parseInt(result.confidence) > 80 ? 'text-warning' : 'text-danger';
        
        return `
            <tr class="result-row" data-result="${index + 1}" style="border-left: 4px solid ${resultColors[index]};">
                <td>
                    <span class="badge" style="background-color: ${resultColors[index]}; color: white;">${result.rank}</span>
                </td>
                <td>
                    <div class="gene-combination-compact">
                        ${result.genes.map(gene => `<span class="badge bg-light text-dark me-1" style="font-size: 0.7rem;">${gene}</span>`).join('')}
                    </div>
                </td>
                <td>
                    <span class="fw-bold ${scoreClass}">${result.score}%</span>
                </td>
                <td>
                    <span class="fw-bold ${confidenceClass}">${result.confidence}%</span>
                </td>
            </tr>
        `;
    }).join('');
    
    tableBody.innerHTML = html;
    
    
    document.getElementById('resultsCard').style.display = 'block';
}

// Toggle heatmap view between grid and collapsed (single row)
function toggleHeatmapView(containerId) {
    let rerenderCallback;
    
    if (containerId === 'optimizeHeatmap') {
        rerenderCallback = renderOptimizeHeatmap;
    }
    
    // Use the common toggle function from common.js
    toggleHeatmapView(containerId, rerenderCallback);
}

