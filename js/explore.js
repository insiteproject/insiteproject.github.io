// Global variables
let datasets = [];
let filteredDatasets = [];
let allUmapData = [];
let umapPlot = null;
let genes = [];
let displayGenes = []; // Limited subset of genes for display
let currentExpression = {};
let currentUncertainty = {};
let filterRowCount = 1;
let maxFilterRows = 2;
let allFilterExpressions = []; // Store expression data for each filter set
let allFilterUncertainties = []; // Store uncertainty data for each filter set

// Performance limits
const MAX_DISPLAY_GENES = 250;  // Max genes to show in heatmap
const MAX_UMAP_POINTS_PER_SAMPLE = 50;  // Max UMAP points per sample
const MAX_TOTAL_UMAP_POINTS = 20000;  // Max total UMAP points

// Load datasets on page load
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
});

// Load datasets using DataService
async function loadDatasets() {
    try {
        datasets = await dataService.loadDatasets();
        genes = await dataService.loadGenes();

        // Keep gene order identical to highly_variable_genes_list.json because
        // mean_expression_profile is stored as an array aligned to that order.
        // Limit genes displayed to MAX_DISPLAY_GENES for performance.
        displayGenes = genes.slice(0, MAX_DISPLAY_GENES);
        console.log('Loaded', datasets.length, 'datasets and', genes.length, 'genes (displaying', displayGenes.length, ')');

        await loadAllUmapData();
        await populateFilters(); // Wait for filters to be populated
        initializeUmapPlot();
        await applyFilters(); // Wait for filters to be applied and expression computed
        // updateVisualization is called inside applyFilters, no need to call again
    } catch (error) {
        console.error('Error loading datasets:', error);
    } finally {
        hideSpinner();
    }
}

// Populate filter dropdowns using DataService
async function populateFilters() {
    try {
        const perturbationTypes = await dataService.getUniqueValues('perturbation_type');
        const tcellSubtypes = await dataService.getUniqueValues('t_cell_subtype');
        const donorTypes = await dataService.getUniqueValues('donor_type');
        const timePoints = await dataService.getUniqueValues('time_point');
        const tissueOrgans = await dataService.getUniqueValues('location');

        // Populate first row (always exists)
        populateSelect('perturbationType_0', perturbationTypes);
        populateSelect('tcellSubtype_0', tcellSubtypes);
        populateSelect('donorType_0', donorTypes);
        populateSelect('timePoint_0', timePoints);
        populateSelect('tissueOrgan_0', tissueOrgans);

        // Add smart filtering - disable unavailable options
        addFilterChangeListeners(0);
        await updateFilterOptions(0);
    } catch (error) {
        console.error('Error populating filters:', error);
    }
}

// Helper function to populate select elements
function populateSelect(elementId, options) {
    const select = document.getElementById(elementId);
    if (!select) return; // Skip if element doesn't exist
    options.forEach(option => {
        const optionElement = document.createElement('option');
        optionElement.value = option;
        optionElement.textContent = option;
        select.appendChild(optionElement);
    });
}

// Update filter options based on current selections - disable unavailable combinations
async function updateFilterOptions(rowIndex = 0) {
    // Get current filter values
    const currentFilters = {
        perturbation_type: document.getElementById(`perturbationType_${rowIndex}`)?.value || 'all',
        t_cell_subtype: document.getElementById(`tcellSubtype_${rowIndex}`)?.value || 'all',
        donor_type: document.getElementById(`donorType_${rowIndex}`)?.value || 'all',
        time_point: document.getElementById(`timePoint_${rowIndex}`)?.value || 'all',
        location: document.getElementById(`tissueOrgan_${rowIndex}`)?.value || 'all'
    };

    // Define filter config: selectId -> dataProperty
    const filterConfig = [
        { selectId: `perturbationType_${rowIndex}`, property: 'perturbation_type' },
        { selectId: `tcellSubtype_${rowIndex}`, property: 't_cell_subtype' },
        { selectId: `donorType_${rowIndex}`, property: 'donor_type' },
        { selectId: `tissueOrgan_${rowIndex}`, property: 'location' }
    ];

    // Update each dropdown
    for (const config of filterConfig) {
        const select = document.getElementById(config.selectId);
        if (!select) continue;

        // Get available values given other current filters
        const availableValues = await dataService.getAvailableValues(config.property, currentFilters);

        // Update option states
        Array.from(select.options).forEach(option => {
            if (option.value === 'all') {
                option.disabled = false;
                option.style.color = '';
            } else {
                const isAvailable = availableValues.has(option.value);
                option.disabled = !isAvailable;
                option.style.color = isAvailable ? '' : '#ccc';
            }
        });
    }
}

// Add change listeners to filter dropdowns for smart filtering
function addFilterChangeListeners(rowIndex = 0) {
    const selectIds = [
        `perturbationType_${rowIndex}`,
        `tcellSubtype_${rowIndex}`,
        `donorType_${rowIndex}`,
        `tissueOrgan_${rowIndex}`
    ];

    selectIds.forEach(selectId => {
        const select = document.getElementById(selectId);
        if (select) {
            select.addEventListener('change', () => updateFilterOptions(rowIndex));
        }
    });
}

// Add a new filter row
function addFilterRow() {
    if (filterRowCount >= maxFilterRows) return;
    
    const filterRows = document.getElementById('filterRows');
    const rowIndex = filterRowCount;
    const color = rowIndex === 1 ? 'blue' : 'green'; // Second row is blue
    
    const newRow = document.createElement('div');
    newRow.className = 'filter-row mb-3';
    newRow.setAttribute('data-row', rowIndex);
    newRow.setAttribute('data-color', color);
    
    newRow.innerHTML = `
        <div class="d-flex flex-wrap gap-3 align-items-end">
            <div class="flex-fill">
                <label class="form-label small">Condition / Perturbation:</label>
                <select class="form-select" id="perturbationType_${rowIndex}">
                    <option value="all">All</option>
                </select>
            </div>
            <div class="flex-fill">
                <label class="form-label small">T-cell:</label>
                <select class="form-select" id="tcellSubtype_${rowIndex}">
                    <option value="all">All</option>
                </select>
            </div>
            <div class="flex-fill">
                <label class="form-label small">Donor:</label>
                <select class="form-select" id="donorType_${rowIndex}">
                    <option value="all">All</option>
                </select>
            </div>
            <div class="flex-fill" style="display: none;">
                <label class="form-label small">Time:</label>
                <select class="form-select" id="timePoint_${rowIndex}">
                    <option value="all">All</option>
                </select>
            </div>
            <div class="flex-fill">
                <label class="form-label small">Tissue:</label>
                <select class="form-select" id="tissueOrgan_${rowIndex}">
                    <option value="all">All</option>
                </select>
            </div>
            <div class="flex-shrink-0" style="width: 80px;">
                <div class="d-flex gap-1">
                    <button class="btn btn-success btn-sm" onclick="addFilterRow()" style="visibility: hidden;" disabled>
                        <i class="fas fa-plus"></i>
                    </button>
                    <button class="btn btn-danger btn-sm" onclick="removeFilterRow(${rowIndex})">
                        <i class="fas fa-minus"></i>
                    </button>
                </div>
            </div>
        </div>
    `;
    
    filterRows.appendChild(newRow);
    filterRowCount++;
    
    // Hide the add button if we've reached the maximum
    if (filterRowCount >= maxFilterRows) {
        document.getElementById('addRowBtn').style.display = 'none';
    }
    
    // Populate the new row's dropdowns
    populateNewRowFilters(rowIndex);
}

// Remove a filter row
function removeFilterRow(rowIndex) {
    const row = document.querySelector(`[data-row="${rowIndex}"]`);
    if (row && rowIndex > 0) { // Don't allow removing the first row
        row.remove();
        filterRowCount--;
        
        // Show the add button again
        document.getElementById('addRowBtn').style.display = 'inline-block';
        
        // Don't auto-apply filters - wait for user to click Apply button
    }
}

// Populate filters for a new row
async function populateNewRowFilters(rowIndex) {
    try {
        const perturbationTypes = await dataService.getUniqueValues('perturbation_type');
        const tcellSubtypes = await dataService.getUniqueValues('t_cell_subtype');
        const donorTypes = await dataService.getUniqueValues('donor_type');
        const timePoints = await dataService.getUniqueValues('time_point');
        const tissueOrgans = await dataService.getUniqueValues('location');

        populateSelect(`perturbationType_${rowIndex}`, perturbationTypes);
        populateSelect(`tcellSubtype_${rowIndex}`, tcellSubtypes);
        populateSelect(`donorType_${rowIndex}`, donorTypes);
        populateSelect(`timePoint_${rowIndex}`, timePoints);
        populateSelect(`tissueOrgan_${rowIndex}`, tissueOrgans);

        // Add smart filtering for the new row
        addFilterChangeListeners(rowIndex);
        await updateFilterOptions(rowIndex);
    } catch (error) {
        console.error('Error populating new row filters:', error);
    }
}

// Apply filters using DataService
async function applyFilters() {
    showSpinner('Applying filters and computing expression profiles...');
    
    try {
        // Collect filters from all active rows
        const allFilterSets = [];
        const filterRows = document.querySelectorAll('.filter-row');
        
        filterRows.forEach((row, index) => {
            const rowIndex = row.getAttribute('data-row');
            const filters = {
                perturbation_type: document.getElementById(`perturbationType_${rowIndex}`)?.value || 'all',
                t_cell_subtype: document.getElementById(`tcellSubtype_${rowIndex}`)?.value || 'all',
                donor_type: document.getElementById(`donorType_${rowIndex}`)?.value || 'all',
                time_point: document.getElementById(`timePoint_${rowIndex}`)?.value || 'all',
                location: document.getElementById(`tissueOrgan_${rowIndex}`)?.value || 'all'
            };
            allFilterSets.push(filters);
        });

        // Use the first filter set for the main filtering logic and summary
        const primaryFilters = allFilterSets[0];
        filteredDatasets = await dataService.filterDatasets(primaryFilters);
        
        // Compute expression data for all filter sets
        await computeAllFilterExpressions(allFilterSets);
        
        updateSelectionSummary();
        updateVisualization(allFilterSets);
        updateDatasetList();
        
    } catch (error) {
        console.error('Error applying filters:', error);
        alert('Error applying filters: ' + error.message);
    } finally {
        hideSpinner();
        
        // Force hide spinner after a short delay as backup
        setTimeout(() => {
            hideSpinner();
        }, 1000);
    }
}

// Update selection summary using DataService
async function updateSelectionSummary() {
    try {
        const stats = await dataService.getStatistics(filteredDatasets);
        document.getElementById('datasetCount').textContent = stats.recordCount;
        document.getElementById('cellCount').textContent = stats.totalCells.toLocaleString();
    } catch (error) {
        console.error('Error updating selection summary:', error);
    }
}

async function loadUnifiedUmapData() {
    try {
        const response = await fetch(`${dataService.baseUrl}data/unified_umap.json`);
        if (!response.ok) {
            return null;
        }
        const payload = await response.json();
        if (Array.isArray(payload)) {
            return payload;
        }
        if (!payload || !Array.isArray(payload.points) || !Array.isArray(payload.fields)) {
            return null;
        }
        const fields = payload.fields;
        const idx = {
            x: fields.indexOf('x'),
            y: fields.indexOf('y'),
            t_cell_subtype: fields.indexOf('t_cell_subtype'),
            donor_type: fields.indexOf('donor_type'),
            location: fields.indexOf('location'),
            perturbation_type: fields.indexOf('perturbation_type'),
            record_name: fields.indexOf('record_name')
        };

        const points = payload.points.map(row => ({
            x: row[idx.x],
            y: row[idx.y],
            t_cell_subtype: row[idx.t_cell_subtype],
            donor_type: row[idx.donor_type],
            location: row[idx.location],
            perturbation_type: row[idx.perturbation_type],
            record_name: row[idx.record_name]
        }));

        if (points.length > MAX_TOTAL_UMAP_POINTS) {
            const step = Math.ceil(points.length / MAX_TOTAL_UMAP_POINTS);
            return points.filter((_, i) => i % step === 0).slice(0, MAX_TOTAL_UMAP_POINTS);
        }
        return points;
    } catch (error) {
        console.warn('Unified UMAP not available:', error);
        return null;
    }
}

// Load all UMAP data from datasets with performance limits
async function loadAllUmapData() {
    allUmapData = [];

    const unified = await loadUnifiedUmapData();
    if (unified && unified.length > 0) {
        allUmapData = unified;
        console.log('Loaded', allUmapData.length, 'unified UMAP points');
        return;
    }

    datasets.forEach(dataset => {
        if (dataset.subsampled_umap && dataset.subsampled_umap.length > 0) {
            // Limit points per sample
            const points = dataset.subsampled_umap;
            const maxPoints = Math.min(points.length, MAX_UMAP_POINTS_PER_SAMPLE);

            // Subsample if needed (take evenly spaced points)
            const step = points.length > maxPoints ? Math.floor(points.length / maxPoints) : 1;

            for (let i = 0; i < points.length && allUmapData.length < MAX_TOTAL_UMAP_POINTS; i += step) {
                const point = points[i];
                allUmapData.push({
                    x: point[0],
                    y: point[1],
                    dataset_id: dataset.id,
                    perturbation_type: dataset.perturbation_type || dataset.perturbation,
                    perturbation: dataset.perturbation,
                    t_cell_subtype: dataset.t_cell_subtype,
                    donor_type: dataset.donor_type,
                    time_point: dataset.time_point,
                    location: dataset.location,
                    record_name: dataset.record_name
                });
            }

            // Stop if we've reached the total limit
            if (allUmapData.length >= MAX_TOTAL_UMAP_POINTS) {
                console.log('Reached UMAP point limit:', MAX_TOTAL_UMAP_POINTS);
                return;
            }
        }
    });

    console.log('Loaded', allUmapData.length, 'UMAP points');
}

// Initialize UMAP plot with all data points in gray
function initializeUmapPlot() {
    if (allUmapData.length === 0) {
        console.warn('No UMAP data available');
        return;
    }

    const trace = {
        x: allUmapData.map(d => d.x),
        y: allUmapData.map(d => d.y),
        mode: 'markers',
        type: 'scatter',
        marker: {
            color: '#d3d3d3', // Pale gray
            size: 3,
            opacity: 0.6
        },
        text: allUmapData.map(d => `${d.record_name || 'Unknown'}<br>Type: ${d.t_cell_subtype}<br>Donor: ${d.donor_type}<br>Location: ${d.location || 'Unknown'}`),
        hovertemplate: '%{text}<extra></extra>',
        showlegend: false
    };

    const layout = {
        title: 'UMAP Visualization',
        xaxis: { title: 'UMAP 1' },
        yaxis: { title: 'UMAP 2' },
        hovermode: 'closest',
        plot_bgcolor: 'white',
        paper_bgcolor: 'white',
        width: 800,
        height: 500,
        dragmode: 'pan'
    };

    const config = {
        responsive: false,
        displayModeBar: true,
        modeBarButtonsToRemove: ['select2d', 'lasso2d', 'autoScale2d'],
        scrollZoom: true
    };

    Plotly.newPlot('umapPlot', [trace], layout, config);
    umapPlot = document.getElementById('umapPlot');
}

// Update UMAP plot to highlight filtered data with multiple colors
function updateUmapPlot(allFilterSets = null) {
    if (!umapPlot || allUmapData.length === 0) {
        return;
    }
    
    // If no filter sets provided, get them from the current rows
    if (!allFilterSets) {
        allFilterSets = [];
        const filterRows = document.querySelectorAll('.filter-row');
        filterRows.forEach((row) => {
            const rowIndex = row.getAttribute('data-row');
            const filters = {
                perturbation_type: document.getElementById(`perturbationType_${rowIndex}`)?.value || 'all',
                t_cell_subtype: document.getElementById(`tcellSubtype_${rowIndex}`)?.value || 'all',
                donor_type: document.getElementById(`donorType_${rowIndex}`)?.value || 'all',
                time_point: document.getElementById(`timePoint_${rowIndex}`)?.value || 'all',
                location: document.getElementById(`tissueOrgan_${rowIndex}`)?.value || 'all'
            };
            allFilterSets.push(filters);
        });
    }
    
    // Define colors for each filter set
    const filterColors = ['#ff0000', '#0000ff']; // Red for first, blue for second
    
    let colors = allUmapData.map(() => '#d3d3d3'); // Default gray
    let sizes = allUmapData.map(() => 3); // Default size
    
    // Check each point against each filter set
    allUmapData.forEach((point, pointIndex) => {
        for (let filterIndex = 0; filterIndex < allFilterSets.length; filterIndex++) {
            const filters = allFilterSets[filterIndex];
            
            // Check if any filters are actually applied (not all set to "all")
            const hasActiveFilters = Object.values(filters).some(value => value !== 'all');
            
            if (hasActiveFilters) {
                const matches = Object.entries(filters).every(([key, value]) => {
                    if (value === 'all') return true;
                    
                    // Map filter keys to actual dataset properties
                    let datasetValue;
                    if (key === 'perturbation_type') {
                        datasetValue = point.perturbation_type || point.perturbation;
                    } else if (key === 'location') {
                        datasetValue = point.location;
                    } else {
                        datasetValue = point[key];
                    }
                    
                    return datasetValue === value;
                });
                
                if (matches) {
                    colors[pointIndex] = filterColors[filterIndex];
                    sizes[pointIndex] = 4;
                    break; // Use the first matching filter's color
                }
            }
        }
    });

    // Update the plot
    const update = {
        'marker.color': [colors],
        'marker.size': sizes
    };

    Plotly.restyle('umapPlot', update, [0]);
}

// Update visualization
function updateVisualization(allFilterSets = null) {
    updateUmapPlot(allFilterSets);
    renderExploreHeatmaps();
    
    // Show/hide and update volcano plot based on number of filter rows
    updateVolcanoVisualization(allFilterSets);
}

// Update volcano visualization based on number of filter rows
async function updateVolcanoVisualization(allFilterSets = null) {
    const volcanoCard = document.getElementById('volcanoCard');
    
    if (!allFilterSets) {
        allFilterSets = [];
        const filterRows = document.querySelectorAll('.filter-row');
        filterRows.forEach((row) => {
            const rowIndex = row.getAttribute('data-row');
            const filters = {
                perturbation_type: document.getElementById(`perturbationType_${rowIndex}`)?.value || 'all',
                t_cell_subtype: document.getElementById(`tcellSubtype_${rowIndex}`)?.value || 'all',
                donor_type: document.getElementById(`donorType_${rowIndex}`)?.value || 'all',
                time_point: document.getElementById(`timePoint_${rowIndex}`)?.value || 'all',
                location: document.getElementById(`tissueOrgan_${rowIndex}`)?.value || 'all'
            };
            allFilterSets.push(filters);
        });
    }
    
    if (allFilterSets.length === 2) {
        // Show volcano plot and perform DE analysis
        volcanoCard.style.display = 'block';
        await performDifferentialExpression(allFilterSets[0], allFilterSets[1]);
    } else {
        // Hide volcano plot
        volcanoCard.style.display = 'none';
    }
}

// Perform differential expression analysis between two filter sets
async function performDifferentialExpression(filterSet1, filterSet2) {
    try {
        // Get datasets for each filter set to access pseudobulk data
        const datasets1 = await dataService.filterDatasets(filterSet1);
        const datasets2 = await dataService.filterDatasets(filterSet2);
        
        // Collect pseudobulk expression data
        const pseudobulk1 = [];
        const pseudobulk2 = [];
        
        datasets1.forEach(dataset => {
            if (dataset.pseudobulk_expression) {
                dataset.pseudobulk_expression.forEach(sample => {
                    pseudobulk1.push(sample.expression);
                });
            }
        });
        
        datasets2.forEach(dataset => {
            if (dataset.pseudobulk_expression) {
                dataset.pseudobulk_expression.forEach(sample => {
                    pseudobulk2.push(sample.expression);
                });
            }
        });
        
        if (pseudobulk1.length === 0 || pseudobulk2.length === 0) {
            console.warn('No pseudobulk data available for DE analysis');
            return;
        }
        
        // Perform DE analysis
        const deResults = performDEAnalysis(pseudobulk1, pseudobulk2, genes);
        
        // Create volcano plot
        createVolcanoPlot(deResults, filterSet1, filterSet2);
        
        // Update description
        updateVolcanoDescription(filterSet1, filterSet2);
        
    } catch (error) {
        console.error('Error performing differential expression analysis:', error);
    }
}

// Perform statistical DE analysis
function performDEAnalysis(group1, group2, geneNames) {
    const results = [];
    
    for (let geneIndex = 0; geneIndex < geneNames.length; geneIndex++) {
        const gene = geneNames[geneIndex];
        
        // Extract expression values for this gene across all samples
        const values1 = group1.map(sample => sample[geneIndex] || 0);
        const values2 = group2.map(sample => sample[geneIndex] || 0);
        
        // Calculate means
        const mean1 = values1.reduce((a, b) => a + b, 0) / values1.length;
        const mean2 = values2.reduce((a, b) => a + b, 0) / values2.length;
        
        // Calculate log2 fold change
        const log2FC = Math.log2((mean2 + 1) / (mean1 + 1));
        
        // Perform t-test
        const pValue = performTTest(values1, values2);
        const negLog10P = -Math.log10(Math.max(pValue, 1e-300)); // Avoid log(0)
        
        results.push({
            gene: gene,
            log2FoldChange: log2FC,
            pValue: pValue,
            negLog10PValue: negLog10P,
            mean1: mean1,
            mean2: mean2
        });
    }
    
    return results;
}

// Proper t-test implementation using stdlib-js
function performTTest(group1, group2) {
    const n1 = group1.length;
    const n2 = group2.length;
    
    if (n1 < 2 || n2 < 2) return 1.0; // Not enough samples
    
    try {
        // Use the stdlib-js ttest function for two-sample t-test
        if (window.ttest) {
            const result = window.ttest(group1, group2);
            return result.pValue;
        } else {
            // Fallback to simple implementation if library not loaded
            return performSimpleTTest(group1, group2);
        }
    } catch (error) {
        console.warn('Error performing t-test, using fallback:', error);
        return performSimpleTTest(group1, group2);
    }
}

// Fallback simple t-test implementation
function performSimpleTTest(group1, group2) {
    const n1 = group1.length;
    const n2 = group2.length;
    
    const mean1 = group1.reduce((a, b) => a + b, 0) / n1;
    const mean2 = group2.reduce((a, b) => a + b, 0) / n2;
    
    const var1 = group1.reduce((sum, x) => sum + Math.pow(x - mean1, 2), 0) / (n1 - 1);
    const var2 = group2.reduce((sum, x) => sum + Math.pow(x - mean2, 2), 0) / (n2 - 1);
    
    const pooledSE = Math.sqrt(var1 / n1 + var2 / n2);
    
    if (pooledSE === 0) return 1.0; // No variance
    
    const tStat = Math.abs(mean1 - mean2) / pooledSE;
    
    // Simple approximation for p-value
    const pValue = Math.min(1.0, 2 * Math.exp(-0.717 * tStat - 0.416 * tStat * tStat));
    
    return pValue;
}

// Create volcano plot using Plotly
function createVolcanoPlot(deResults, filterSet1, filterSet2) {
    const significanceThreshold = 0.05;
    const foldChangeThreshold = 1.0; // log2 fold change threshold
    
    // Separate points by significance and fold change
    const significant = deResults.filter(d => d.pValue < significanceThreshold && Math.abs(d.log2FoldChange) > foldChangeThreshold);
    const nonSignificant = deResults.filter(d => d.pValue >= significanceThreshold || Math.abs(d.log2FoldChange) <= foldChangeThreshold);
    
    const traces = [
        {
            x: nonSignificant.map(d => d.log2FoldChange),
            y: nonSignificant.map(d => d.negLog10PValue),
            mode: 'markers',
            type: 'scatter',
            name: 'Non-significant',
            marker: {
                color: '#d3d3d3',
                size: 4,
                opacity: 0.6
            },
            text: nonSignificant.map(d => `${d.gene}<br>Log2FC: ${d.log2FoldChange.toFixed(2)}<br>P-value: ${d.pValue.toExponential(2)}`),
            hovertemplate: '%{text}<extra></extra>',
            showlegend: true
        },
        {
            x: significant.map(d => d.log2FoldChange),
            y: significant.map(d => d.negLog10PValue),
            mode: 'markers',
            type: 'scatter',
            name: 'Significant',
            marker: {
                color: '#ff0000',
                size: 6,
                opacity: 0.8
            },
            text: significant.map(d => `${d.gene}<br>Log2FC: ${d.log2FoldChange.toFixed(2)}<br>P-value: ${d.pValue.toExponential(2)}`),
            hovertemplate: '%{text}<extra></extra>',
            showlegend: true
        }
    ];
    
    const layout = {
        title: 'Volcano Plot - Differential Expression',
        xaxis: { 
            title: 'Log2 Fold Change (Filter 2 vs Filter 1)',
            zeroline: true,
            zerolinecolor: '#000000',
            zerolinewidth: 1
        },
        yaxis: { 
            title: '-Log10 P-value',
            zeroline: false
        },
        hovermode: 'closest',
        plot_bgcolor: 'white',
        paper_bgcolor: 'white',
        shapes: [
            // Horizontal line for significance threshold
            {
                type: 'line',
                x0: -10,
                x1: 10,
                y0: -Math.log10(significanceThreshold),
                y1: -Math.log10(significanceThreshold),
                line: {
                    color: 'red',
                    width: 1,
                    dash: 'dash'
                }
            },
            // Vertical lines for fold change thresholds
            {
                type: 'line',
                x0: foldChangeThreshold,
                x1: foldChangeThreshold,
                y0: 0,
                y1: 20,
                line: {
                    color: 'red',
                    width: 1,
                    dash: 'dash'
                }
            },
            {
                type: 'line',
                x0: -foldChangeThreshold,
                x1: -foldChangeThreshold,
                y0: 0,
                y1: 20,
                line: {
                    color: 'red',
                    width: 1,
                    dash: 'dash'
                }
            }
        ]
    };
    
    const config = {
        responsive: false,
        displayModeBar: true,
        modeBarButtonsToRemove: ['select2d', 'lasso2d', 'autoScale2d'],
        scrollZoom: true
    };
    
    Plotly.newPlot('volcanoPlotContainer', traces, layout, config);
}

// Update volcano plot description with filter details
function updateVolcanoDescription(filterSet1, filterSet2) {
    const descriptionElement = document.getElementById('volcanoDescription');
    if (!descriptionElement) return;
    
    // Helper function to format filter set for display
    function formatFilterSet(filters) {
        const activeFilters = [];
        Object.entries(filters).forEach(([key, value]) => {
            if (value !== 'all') {
                const displayKey = key.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase());
                activeFilters.push(`${displayKey}: ${value}`);
            }
        });
        return activeFilters.length > 0 ? activeFilters.join(', ') : 'All conditions';
    }
    
    const filter1Description = formatFilterSet(filterSet1);
    const filter2Description = formatFilterSet(filterSet2);
    
    descriptionElement.innerHTML = `
        Differential gene expression analysis comparing two filter sets:<br>
        <strong>Filter 1 (Red):</strong> ${filter1Description}<br>
        <strong>Filter 2 (Blue):</strong> ${filter2Description}<br>
        Points represent genes, with log2 fold change (Filter 2 vs Filter 1) on x-axis and -log10 p-value on y-axis.
    `;
}


// Update dataset list
function updateDatasetList() {
    const datasetList = document.getElementById('datasetList');
    
    if (filteredDatasets.length === 0) {
        datasetList.innerHTML = '<p class="text-muted">No records match the current filters.</p>';
        return;
    }
    
    // Group datasets by parent record to avoid duplicates
    const parentRecords = new Map();
    
    filteredDatasets.forEach(dataset => {
        const recordKey = `${dataset.repo_name}_${dataset.id_on_repo}`;
        if (!parentRecords.has(recordKey)) {
            parentRecords.set(recordKey, {
                record_name: dataset.record_name,
                repo_name: dataset.repo_name,
                id_on_repo: dataset.id_on_repo,
                url: dataset.url,
                description: dataset.description,
                totalDataObjects: 0,
                totalCells: 0
            });
        }
        
        const record = parentRecords.get(recordKey);
        record.totalDataObjects++;
        record.totalCells += dataset.cell_count || 0;
    });
    
    const html = Array.from(parentRecords.values()).map(record => `
        <div class="card mb-2">
            <div class="card-body">
                <h6 class="card-title">${record.record_name}</h6>
                <div class="row">
                    <div class="col-md-3"><small><strong>Repository:</strong> ${record.repo_name}</small></div>
                    <div class="col-md-3"><small><strong>Data Objects:</strong> ${record.totalDataObjects}</small></div>
                    <div class="col-md-3"><small><strong>Total Cells:</strong> ${record.totalCells.toLocaleString()}</small></div>
                    <div class="col-md-3"><small><strong>URL:</strong> <a href="${record.url}" target="_blank" class="text-decoration-none">View Record</a></small></div>
                </div>
            </div>
        </div>
    `).join('');
    
    datasetList.innerHTML = html;
}

// Compute expression values from filtered datasets
async function computeCurrentExpression(filters) {
    try {
        const expressionProfile = await dataService.computeExpressionProfile(filters);
        
        genes.forEach(gene => {
            if (expressionProfile[gene]) {
                currentExpression[gene] = expressionProfile[gene].mean;
                currentUncertainty[gene] = expressionProfile[gene].std;
            } else {
                // Fallback to zero if gene not found
                currentExpression[gene] = 0;
                currentUncertainty[gene] = 0;
            }
        });
        
    } catch (error) {
        console.error('Error computing current expression:', error);
        // Fallback to zero values
        genes.forEach(gene => {
            currentExpression[gene] = 0;
            currentUncertainty[gene] = 0;
        });
    }
}

// Compute expression values for all filter sets
async function computeAllFilterExpressions(allFilterSets) {
    allFilterExpressions = [];
    allFilterUncertainties = [];

    console.log('computeAllFilterExpressions called with', allFilterSets.length, 'filter sets');
    console.log('Global genes array length:', genes.length);

    for (let i = 0; i < allFilterSets.length; i++) {
        const filters = allFilterSets[i];

        try {
            const expressionProfile = await dataService.computeExpressionProfile(filters);

            // Debug: check the expression profile
            const profileGenes = Object.keys(expressionProfile);
            console.log(`Filter ${i}: expressionProfile has ${profileGenes.length} genes`);
            if (profileGenes.length > 0) {
                const sampleGene = profileGenes[0];
                console.log(`Sample: ${sampleGene} = `, expressionProfile[sampleGene]);
            }

            const expression = {};
            const uncertainty = {};

            genes.forEach(gene => {
                if (expressionProfile[gene]) {
                    expression[gene] = expressionProfile[gene].mean;
                    uncertainty[gene] = expressionProfile[gene].std;
                } else {
                    expression[gene] = 0;
                    uncertainty[gene] = 0;
                }
            });

            // Debug: check computed values
            if (genes.length > 0) {
                console.log(`Expression for ${genes[0]}: ${expression[genes[0]]}`);
            }

            allFilterExpressions.push(expression);
            allFilterUncertainties.push(uncertainty);
            
        } catch (error) {
            console.error(`Error computing expression for filter set ${i}:`, error);
            // Fallback to zero values
            const expression = {};
            const uncertainty = {};
            genes.forEach(gene => {
                expression[gene] = 0;
                uncertainty[gene] = 0;
            });
            allFilterExpressions.push(expression);
            allFilterUncertainties.push(uncertainty);
        }
    }
    
    // Update the primary expression data with the first filter set
    if (allFilterExpressions.length > 0) {
        currentExpression = allFilterExpressions[0];
        currentUncertainty = allFilterUncertainties[0];
    }
}

// Render expression heatmaps for explore page
function renderExploreHeatmaps() {
    const container = document.getElementById('expressionProfilesContainer');
    if (!container || displayGenes.length === 0) return;
    
    // Clear existing content
    container.innerHTML = '';
    
    // Define colors and labels for each filter set
    const filterColors = ['red', 'blue'];
    const filterLabels = ['Filter 1', 'Filter 2'];
    
    // Create separate cards for each active filter set
    for (let i = 0; i < allFilterExpressions.length; i++) {
        const expression = allFilterExpressions[i];
        const uncertainty = allFilterUncertainties[i];
        const color = filterColors[i] || 'gray';
        const label = filterLabels[i] || `Filter ${i + 1}`;
        
        // Create card for this heatmap
        const cardDiv = document.createElement('div');
        cardDiv.className = 'card card-surface mb-4';
        const geneCountInfo = displayGenes.length < genes.length ?
            `<small class="text-muted ms-2">(showing ${displayGenes.length} of ${genes.length} genes)</small>` : '';
        cardDiv.innerHTML = `
            <div class="card-header d-flex justify-content-between align-items-center">
                <h5 style="color: ${color};">
                    <i class="fas fa-circle" style="color: ${color}; font-size: 0.8em;"></i>
                    ${label} Expression Profile ${geneCountInfo}
                </h5>

            </div>
            <div class="card-body">
                <p class="text-muted mb-3">
                    <small id="filterDescription_${i}">Gene expression profile for the selected cell populations. Expression values are weighted averages across all matching datasets.</small>
                </p>
                <div id="exploreHeatmap_${i}" class="heatmap-container"></div>
            </div>
        `;
        
        container.appendChild(cardDiv);
        
        // Render the heatmap (using displayGenes for performance)
        const heatmapContainer = document.getElementById(`exploreHeatmap_${i}`);
        renderHeatmap(heatmapContainer, displayGenes, expression, uncertainty, {
            interactive: false
        });
        
        // Update filter description
        updateFilterDescription(i);
    }
}

// Update filter description for a specific filter index
function updateFilterDescription(filterIndex) {
    const descriptionElement = document.getElementById(`filterDescription_${filterIndex}`);
    if (!descriptionElement) return;
    
    const filterRows = document.querySelectorAll('.filter-row');
    if (filterIndex >= filterRows.length) return;
    
    const row = filterRows[filterIndex];
    const rowIndex = row.getAttribute('data-row');
    
    const filters = {
        perturbation_type: document.getElementById(`perturbationType_${rowIndex}`)?.value || 'all',
        t_cell_subtype: document.getElementById(`tcellSubtype_${rowIndex}`)?.value || 'all',
        donor_type: document.getElementById(`donorType_${rowIndex}`)?.value || 'all',
        time_point: document.getElementById(`timePoint_${rowIndex}`)?.value || 'all',
        location: document.getElementById(`tissueOrgan_${rowIndex}`)?.value || 'all'
    };
    
    // Format filter description
    const activeFilters = [];
    Object.entries(filters).forEach(([key, value]) => {
        if (value !== 'all') {
            const displayKey =
                key === 'perturbation_type' ? 'Condition' :
                key.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase());
            activeFilters.push(`${displayKey}: ${value}`);
        }
    });
    
    const description = activeFilters.length > 0 ? 
        `Filters: ${activeFilters.join(', ')}. Gene expression profile for the selected cell populations. Expression values are weighted averages across all matching datasets.` : 
        'All conditions selected. Gene expression profile for all cell populations. Expression values are weighted averages across all datasets.';
    descriptionElement.textContent = description;
}

// Toggle all heatmap views between grid and collapsed (legacy function, now handled per card)
function toggleAllHeatmapViews() {
    const containers = document.querySelectorAll('[id^="exploreHeatmap_"]');
    
    if (containers.length === 0) return;
    
    // Check current state from first container
    const isCollapsed = containers[0].classList.contains('heatmap-collapsed');
    
    containers.forEach((container, index) => {
        const button = document.getElementById(`exploreHeatmapCollapseBtn_${index}`);
        
        if (isCollapsed) {
            // Expand to grid view
            container.classList.remove('heatmap-collapsed');
            if (button) button.innerHTML = '<i class="fas fa-compress-alt"></i> Collapse';
        } else {
            // Collapse to single row
            container.classList.add('heatmap-collapsed');
            if (button) button.innerHTML = '<i class="fas fa-expand-alt"></i> Expand';
        }
    });
    
    // Re-render all heatmaps
    renderExploreHeatmaps();
}
