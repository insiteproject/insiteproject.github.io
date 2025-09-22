// Global variables
let datasets = [];
let filteredDatasets = [];
let allUmapData = [];
let umapPlot = null;
let genes = [];
let currentExpression = {};
let currentUncertainty = {};

// Load datasets on page load
document.addEventListener('DOMContentLoaded', function() {
    loadDatasets();
});

// Load datasets using DataService
async function loadDatasets() {
    try {
        datasets = await dataService.loadDatasets();
        genes = await dataService.loadGenes();
        genes.sort(); // Sort alphabetically
        await loadAllUmapData();
        populateFilters();
        updateVolcanoOptions(); // Initialize volcano options
        initializeUmapPlot();
        applyFilters();
        // Load initial visualization
        updateVisualization();
    } catch (error) {
        console.error('Error loading datasets:', error);
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

        populateSelect('perturbationType', perturbationTypes);
        populateSelect('tcellSubtype', tcellSubtypes);
        populateSelect('donorType', donorTypes);
        populateSelect('timePoint', timePoints);
        populateSelect('tissueOrgan', tissueOrgans);
    } catch (error) {
        console.error('Error populating filters:', error);
    }
}

// Helper function to populate select elements
function populateSelect(elementId, options) {
    const select = document.getElementById(elementId);
    options.forEach(option => {
        const optionElement = document.createElement('option');
        optionElement.value = option;
        optionElement.textContent = option;
        select.appendChild(optionElement);
    });
}

// Apply filters using DataService
async function applyFilters() {
    try {
        const filters = {
            perturbation_type: document.getElementById('perturbationType').value,
            t_cell_subtype: document.getElementById('tcellSubtype').value,
            donor_type: document.getElementById('donorType').value,
            time_point: document.getElementById('timePoint').value,
            location: document.getElementById('tissueOrgan').value
        };

        filteredDatasets = await dataService.filterDatasets(filters);
        await computeCurrentExpression(filters);
        updateSelectionSummary();
        updateVisualization();
        updateDatasetList();
    } catch (error) {
        console.error('Error applying filters:', error);
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

// Load all UMAP data from datasets
async function loadAllUmapData() {
    allUmapData = [];
    
    datasets.forEach(dataset => {
        if (dataset.subsampled_umap && dataset.subsampled_umap.length > 0) {
            dataset.subsampled_umap.forEach(point => {
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
            });
        }
    });
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
        paper_bgcolor: 'white'
    };

    const config = {
        responsive: true,
        displayModeBar: true,
        modeBarButtonsToRemove: ['pan2d', 'select2d', 'lasso2d', 'autoScale2d']
    };

    Plotly.newPlot('umapPlot', [trace], layout, config);
    umapPlot = document.getElementById('umapPlot');
}

// Update UMAP plot to highlight filtered data
function updateUmapPlot() {
    if (!umapPlot || allUmapData.length === 0) {
        return;
    }
    
    // Check if any filters are actually applied (not all set to "all")
    const filters = {
        perturbation_type: document.getElementById('perturbationType').value,
        t_cell_subtype: document.getElementById('tcellSubtype').value,
        donor_type: document.getElementById('donorType').value,
        time_point: document.getElementById('timePoint').value,
        location: document.getElementById('tissueOrgan').value
    };
    
    const hasActiveFilters = Object.values(filters).some(value => value !== 'all');
    
    let colors, sizes;
    
    if (!hasActiveFilters) {
        // No filters applied - show all points in gray
        colors = allUmapData.map(() => '#d3d3d3');
        sizes = allUmapData.map(() => 3);
    } else {
        // Filters applied - highlight matching points by checking each point's properties directly
        colors = allUmapData.map((point) => {
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
            return matches ? '#ffd700' : '#d3d3d3'; // Yellow for filtered, gray for others
        });
        
        sizes = allUmapData.map((point) => {
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
            return matches ? 4 : 3;
        });
    }

    // Update the plot
    const update = {
        'marker.color': [colors],
        'marker.size': sizes
    };

    Plotly.restyle('umapPlot', update, [0]);
}

// Update visualization
function updateVisualization() {
    updateUmapPlot();
    renderExploreHeatmap();
    
    // Also update volcano plot
    updateVolcanoPlot();
}

// Update volcano plot options based on selected variable
async function updateVolcanoOptions() {
    const variable = document.getElementById('volcanoVariable').value;
    const baselineSelect = document.getElementById('volcanoBaseline');
    const targetSelect = document.getElementById('volcanoTarget');
    
    // Clear existing options
    baselineSelect.innerHTML = '';
    targetSelect.innerHTML = '';
    
    let options = [];
    
    // Get options based on selected variable using DataService
    if (variable === 'control_vs_treated') {
        options = ['Control', 'Treated'];
    } else if (variable === 'tcell_subtype') {
        options = await dataService.getUniqueValues('t_cell_subtype');
    } else if (variable === 'time_point') {
        options = await dataService.getUniqueValues('time_point');
    } else if (variable === 'donor_type') {
        options = await dataService.getUniqueValues('donor_type');
    } else if (variable === 'location') {
        options = await dataService.getUniqueValues('location');
    }
    
    // Populate both dropdowns with the same options
    options.forEach(option => {
        const baselineOption = document.createElement('option');
        baselineOption.value = option.toLowerCase().replace(/\s+/g, '_').replace(/\+/g, '');
        baselineOption.textContent = option;
        baselineSelect.appendChild(baselineOption);
        
        const targetOption = document.createElement('option');
        targetOption.value = option.toLowerCase().replace(/\s+/g, '_').replace(/\+/g, '');
        targetOption.textContent = option;
        targetSelect.appendChild(targetOption);
    });
    
    // Set different default values if possible
    if (options.length > 1) {
        targetSelect.selectedIndex = 1;
    }
    
    // Update the plot
    updateVolcanoPlot();
}

// Update volcano plot
function updateVolcanoPlot() {
    const variable = document.getElementById('volcanoVariable').value;
    const baseline = document.getElementById('volcanoBaseline').value;
    const target = document.getElementById('volcanoTarget').value;
    const filters = {
        perturbationType: document.getElementById('perturbationType').value,
        tcellSubtype: document.getElementById('tcellSubtype').value,
        donorType: document.getElementById('donorType').value,
        timePoint: document.getElementById('timePoint').value,
        tissueOrgan: document.getElementById('tissueOrgan').value
    };

    // Load volcano plot with SHA256-based filename
    loadVolcanoImage(filters, variable, baseline, target);
}

// Recursively sort dictionary by keys (matches Python implementation)
function sortDict(obj) {
    if (typeof obj === 'object' && obj !== null && !Array.isArray(obj)) {
        const sorted = {};
        Object.keys(obj).sort().forEach(key => {
            sorted[key] = sortDict(obj[key]);
        });
        return sorted;
    }
    return obj;
}

// Generate SHA256 hash from parameters
async function generateSHA256Hash(params) {
    // Recursively sort all nested objects to match Python implementation
    const sortedParams = sortDict(params);
    
    const jsonString = JSON.stringify(sortedParams);
    
    const encoder = new TextEncoder();
    const data = encoder.encode(jsonString);
    const hashBuffer = await crypto.subtle.digest('SHA-256', data);
    const hashArray = Array.from(new Uint8Array(hashBuffer));
    const hashHex = hashArray.map(b => b.toString(16).padStart(2, '0')).join('');
    return hashHex;
}

// Load volcano image with SHA256-based filename
async function loadVolcanoImage(filters, variable, baseline, target) {
    const params = {
        baseline: baseline,
        filters: filters,
        target: target,
        type: 'volcano',
        variable: variable
    };
    
    const hash = await generateSHA256Hash(params);
    
    // Create query string manually to handle object serialization properly
    const queryParams = new URLSearchParams();
    queryParams.append('filters', JSON.stringify(filters));
    queryParams.append('variable', variable);
    queryParams.append('baseline', baseline);
    queryParams.append('target', target);
    queryParams.append('type', 'volcano');
    
    const imagePath = `images/explore/volcanos/${hash}.png?${queryParams.toString()}`;
    
    document.getElementById('volcanoPlot').src = imagePath;
}

// Update dataset list
function updateDatasetList() {
    const datasetList = document.getElementById('datasetList');
    
    if (filteredDatasets.length === 0) {
        datasetList.innerHTML = '<p class="text-muted">No records match the current filters.</p>';
        return;
    }
    
    const html = filteredDatasets.map(dataset => `
        <div class="card mb-2">
            <div class="card-body">
                <h6 class="card-title">${dataset.name}</h6>
                <div class="row">
                    <div class="col-md-3"><small><strong>Perturbation:</strong> ${dataset.perturbation_type}</small></div>
                    <div class="col-md-3"><small><strong>T-cell type:</strong> ${dataset.t_cell_subtype}</small></div>
                    <div class="col-md-3"><small><strong>Donor:</strong> ${dataset.donor_type}</small></div>
                    <div class="col-md-3"><small><strong>Cells:</strong> ${dataset.cell_count.toLocaleString()}</small></div>
                </div>
                <div class="row mt-1">
                    <div class="col-md-3"><small><strong>Time:</strong> ${dataset.time_point}</small></div>
                    <div class="col-md-3"><small><strong>Tissue:</strong> ${dataset.location}</small></div>
                    <div class="col-md-6"><small><strong>Description:</strong> ${dataset.description}</small></div>
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

// Render expression heatmap for explore page
function renderExploreHeatmap() {
    const container = document.getElementById('exploreHeatmap');
    if (!container || genes.length === 0) return;
    
    // Check if container is in collapsed state
    const isCollapsed = container.classList.contains('heatmap-collapsed');
    const gridClass = isCollapsed ? 'heatmap-grid-collapsed' : 'heatmap-grid';
    
    let html = generateColorbarHTML() + `<div class="${gridClass}">`;
    
    genes.forEach(gene => {
        const expression = currentExpression[gene] || 0;
        const uncertainty = currentUncertainty[gene] || 0;
        const expressionClass = getExpressionClass(expression);
        const uncertaintyClass = getUncertaintyClass(uncertainty);
        
        html += `
            <div class="heatmap-cell ${expressionClass}" 
                 data-gene="${gene}"
                 title="${gene}: Mean=${expression.toFixed(2)}, Std=${uncertainty.toFixed(2)}">
                <div class="heatmap-gene-name">${gene}</div>
                <div class="heatmap-expression">${expression.toFixed(1)}</div>
                <div class="uncertainty-bar ${uncertaintyClass}">
                    <span class="uncertainty-value">${uncertainty.toFixed(1)}</span>
                </div>
            </div>
        `;
    });
    
    html += '</div>';
    container.innerHTML = html;
}

// Get CSS class based on expression level (reused from predict.js)
function getExpressionClass(expression) {
    if (expression < 2) return 'expr-low';
    if (expression < 4) return 'expr-low-med';
    if (expression < 6) return 'expr-medium';
    if (expression < 8) return 'expr-med-high';
    return 'expr-high';
}

// Get CSS class based on uncertainty level (reused from predict.js)
function getUncertaintyClass(uncertainty) {
    if (uncertainty < 0.4) return 'uncert-very-low';
    if (uncertainty < 0.8) return 'uncert-low';
    if (uncertainty < 1.2) return 'uncert-medium';
    if (uncertainty < 1.6) return 'uncert-high';
    return 'uncert-very-high';
}

// Toggle heatmap view between grid and collapsed (reused from predict.js)
function toggleHeatmapView(containerId) {
    const container = document.getElementById(containerId);
    let button;
    
    // Handle different button ID patterns
    if (containerId === 'exploreHeatmap') {
        button = document.getElementById('exploreHeatmapCollapseBtn');
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
    
    // Re-render the heatmap
    renderExploreHeatmap();
}
