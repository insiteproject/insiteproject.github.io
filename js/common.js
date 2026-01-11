// Common JavaScript functionality for INSITE website

/**
 * Data service for managing dataset and gene data access
 */
class DataService {
    constructor() {
        this.datasetsCache = null;
        this.genesCache = null;
        this.baseUrl = '';
    }

    /**
     * Load datasets from the data source
     * @returns {Promise<Array>} Array of dataset objects (flattened data objects)
     */
    async loadDatasets() {
        if (this.datasetsCache) {
            return this.datasetsCache;
        }
        
        try {
            const response = await fetch(`${this.baseUrl}data/datasets.json`);
            if (!response.ok) {
                throw new Error(`HTTP error! status: ${response.status}`);
            }
            const repositories = await response.json();
            
            // Check if this is the new format (array of repositories) or old format (direct data objects)
            if (Array.isArray(repositories) && repositories.length > 0 && repositories[0].data_objects) {
                // New format: flatten data objects from all repositories
                this.datasetsCache = [];
                repositories.forEach(repo => {
                    repo.data_objects.forEach(dataObj => {
                        // Add repository metadata to each data object
                        this.datasetsCache.push({
                            ...dataObj,
                            repo_id: repo.id,
                            repo_name: repo.repo_name,
                            repo_url: repo.url,
                            repo_description: repo.description,
                            record_name: repo.record_name
                        });
                    });
                });
            } else {
                // New direct format: data objects are directly in the array
                this.datasetsCache = repositories;
            }
            
            return this.datasetsCache;
        } catch (error) {
            console.error('Error loading datasets:', error);
            return [];
        }
    }

    /**
     * Load genes from the data source
     * @returns {Promise<Array>} Array of gene names
     */
    async loadGenes() {
        if (this.genesCache) {
            return this.genesCache;
        }
        
        try {
            const response = await fetch(`${this.baseUrl}data/highly_variable_genes_list.json`);
            if (!response.ok) {
                throw new Error(`HTTP error! status: ${response.status}`);
            }
            this.genesCache = await response.json();
            return this.genesCache;
        } catch (error) {
            console.error('Error loading genes:', error);
            return [];
        }
    }

    /**
     * Get unique values for a specific dataset property
     * @param {string} property - Property name (e.g., 't_cell_subtype', 'donor_type')
     * @returns {Promise<Array>} Array of unique values
     */
    async getUniqueValues(property) {
        const datasets = await this.loadDatasets();
        return [...new Set(datasets.map(d => d[property]))].filter(Boolean);
    }

    /**
     * Filter datasets based on criteria
     * @param {Object} filters - Filter criteria object
     * @returns {Promise<Array>} Filtered datasets (individual data objects)
     */
    async filterDatasets(filters) {
        const datasets = await this.loadDatasets();
        return datasets.filter(dataset => {
            return Object.entries(filters).every(([key, value]) => {
                if (value === 'all' || value === null || value === undefined) {
                    return true;
                }
                // Handle the tissue/location field mapping
                if (key === 'location' && dataset.tissue_or_organ) {
                    return dataset.tissue_or_organ === value;
                }
                // Handle perturbation field mapping
                if (key === 'perturbation_type' && dataset.perturbation) {
                    return dataset.perturbation === value;
                }
                return dataset[key] === value;
            });
        });
    }

    /**
     * Get dataset statistics
     * @param {Array} datasets - Optional filtered datasets, uses all if not provided
     * @returns {Promise<Object>} Statistics object
     */
    async getStatistics(datasets = null) {
        if (!datasets) {
            datasets = await this.loadDatasets();
        }

        const totalCells = datasets.reduce((sum, dataset) => sum + (dataset.cell_count || 0), 0);
        const uniquePerturbations = new Set(datasets.map(d => d.perturbation_type)).size;
        const uniqueTCellTypes = new Set(datasets.map(d => d.t_cell_subtype)).size;

        // Count number of data objects (samples/datasets), not unique repositories
        const recordCount = datasets.length;

        return {
            recordCount,
            totalCells,
            uniquePerturbations,
            uniqueTCellTypes
        };
    }

    /**
     * Clear cached data (useful for testing or data refresh)
     */
    clearCache() {
        this.datasetsCache = null;
        this.genesCache = null;
    }

    /**
     * Get available values for a property given current filter selections
     * This enables smart filtering - disabling options that would result in empty sets
     * @param {string} property - The property to get available values for
     * @param {Object} currentFilters - Current filter selections (excluding the target property)
     * @returns {Promise<Set>} Set of available values
     */
    async getAvailableValues(property, currentFilters) {
        const datasets = await this.loadDatasets();

        // Filter datasets based on current selections (excluding the target property)
        const filtered = datasets.filter(dataset => {
            return Object.entries(currentFilters).every(([key, value]) => {
                if (key === property) return true; // Skip the property we're checking
                if (value === 'all' || value === null || value === undefined) return true;
                if (key === 'location' && dataset.tissue_or_organ) {
                    return dataset.tissue_or_organ === value;
                }
                if (key === 'perturbation_type' && dataset.perturbation) {
                    return dataset.perturbation === value;
                }
                return dataset[key] === value;
            });
        });

        // Return unique values of the target property from filtered results
        return new Set(filtered.map(d => d[property]).filter(Boolean));
    }

    /**
     * Set base URL for data requests (useful for different environments)
     * @param {string} url - Base URL
     */
    setBaseUrl(url) {
        this.baseUrl = url.endsWith('/') ? url : url + '/';
    }

    /**
     * Compute weighted mean expression profile from filtered datasets
     * @param {Object} filters - Filter criteria object
     * @returns {Promise<Object>} Object with mean and std arrays for each gene
     */
    async computeExpressionProfile(filters) {
        console.log('Computing expression profile with filters:', filters);
        
        const filteredDatasets = await this.filterDatasets(filters);
        const genes = await this.loadGenes();
        
        console.log('Filtered datasets:', filteredDatasets.length, 'genes:', genes.length);
        
        if (filteredDatasets.length === 0) {
            console.log('No datasets match filters, returning zero expression');
            // Return zero expression if no datasets match
            const result = {};
            genes.forEach(gene => {
                result[gene] = { mean: 0, std: 0 };
            });
            return result;
        }
        
        // Initialize weighted sums
        const weightedMeans = {};
        const weightedStds = {};
        let totalWeight = 0;
        
        genes.forEach(gene => {
            weightedMeans[gene] = 0;
            weightedStds[gene] = 0;
        });
        
        // Compute weighted average across datasets
        let datasetsWithExpression = 0;
        filteredDatasets.forEach((dataset, datasetIndex) => {
            const weight = dataset.n_obs || dataset.cell_count || 1; // Use cell count as weight
            totalWeight += weight;

            // Debug: log first dataset's structure
            if (datasetIndex === 0) {
                console.log('First dataset keys:', Object.keys(dataset));
                console.log('First dataset has mean_expression_profile:', 'mean_expression_profile' in dataset);
                if (dataset.mean_expression_profile) {
                    console.log('Profile length:', dataset.mean_expression_profile.length, 'Expected:', genes.length);
                    console.log('First 3 profile values:', dataset.mean_expression_profile.slice(0, 3));
                }
            }

            if (dataset.mean_expression_profile && dataset.mean_expression_profile.length === genes.length) {
                datasetsWithExpression++;
                genes.forEach((gene, index) => {
                    const profile = dataset.mean_expression_profile[index];
                    if (profile && Array.isArray(profile) && profile.length === 2) {
                        weightedMeans[gene] += profile[0] * weight; // mean is first element
                        weightedStds[gene] += profile[1] * weight;  // std is second element
                    }
                });
            }
        });

        console.log('Datasets with expression data:', datasetsWithExpression, 'total weight:', totalWeight);
        // Debug: log sample results
        if (genes.length > 0) {
            const sampleGene = genes[0];
            console.log(`Sample result for ${sampleGene}: mean=${weightedMeans[sampleGene] / totalWeight}, std=${weightedStds[sampleGene] / totalWeight}`);
        }
        
        // Normalize by total weight
        const result = {};
        genes.forEach(gene => {
            result[gene] = {
                mean: totalWeight > 0 ? weightedMeans[gene] / totalWeight : 0,
                std: totalWeight > 0 ? weightedStds[gene] / totalWeight : 0
            };
        });
        
        console.log('Expression profile computed successfully');
        return result;
    }
}

/**
 * Generate colorbar HTML for heatmaps
 * Thresholds calibrated for log-normalized scRNA-seq expression data
 * @returns {string} HTML string for colorbar
 */
function generateColorbarHTML() {
    return `
        <div class="heatmap-colorbars">
            <div class="colorbar-container">
                <div class="colorbar-label">Expression (Mean)</div>
                <div class="colorbar expression-colorbar">
                    <div class="colorbar-item expr-low">&lt;0.25</div>
                    <div class="colorbar-item expr-low-med">0.25-0.5</div>
                    <div class="colorbar-item expr-medium">0.5-1</div>
                    <div class="colorbar-item expr-med-high">1-2</div>
                    <div class="colorbar-item expr-high">&gt;2</div>
                </div>
            </div>
            <div class="colorbar-container">
                <div class="colorbar-label">Uncertainty (Std)</div>
                <div class="colorbar uncertainty-colorbar">
                    <div class="colorbar-item uncert-very-low">&lt;0.2</div>
                    <div class="colorbar-item uncert-low">0.2-0.4</div>
                    <div class="colorbar-item uncert-medium">0.4-0.7</div>
                    <div class="colorbar-item uncert-high">0.7-1</div>
                    <div class="colorbar-item uncert-very-high">&gt;1</div>
                </div>
            </div>
        </div>`;
}

/**
 * Get CSS class based on expression level
 * Thresholds calibrated for log-normalized scRNA-seq data
 * @param {number} expression - Expression value (log-normalized)
 * @returns {string} CSS class name
 */
function getExpressionClass(expression) {
    if (expression < 0.25) return 'expr-low';
    if (expression < 0.5) return 'expr-low-med';
    if (expression < 1) return 'expr-medium';
    if (expression < 2) return 'expr-med-high';
    return 'expr-high';
}

/**
 * Get CSS class based on uncertainty level
 * Thresholds calibrated for scRNA-seq standard deviation values
 * @param {number} uncertainty - Uncertainty value (std)
 * @returns {string} CSS class name
 */
function getUncertaintyClass(uncertainty) {
    if (uncertainty < 0.2) return 'uncert-very-low';
    if (uncertainty < 0.4) return 'uncert-low';
    if (uncertainty < 0.7) return 'uncert-medium';
    if (uncertainty < 1) return 'uncert-high';
    return 'uncert-very-high';
}

/**
 * Unified heatmap renderer for gene expression data
 * @param {HTMLElement} container - Container element to render heatmap in
 * @param {Array} genes - Array of gene names
 * @param {Object} expressionData - Object mapping gene names to expression values
 * @param {Object} uncertaintyData - Object mapping gene names to uncertainty values
 * @param {Object} options - Rendering options
 * @param {boolean} options.interactive - Whether cells should be clickable
 * @param {Set} options.selectedGenes - Set of selected gene names (for knockout visualization)
 * @param {Function} options.onGeneClick - Callback function for gene clicks
 */
function renderHeatmap(container, genes, expressionData, uncertaintyData, options = {}) {
    const {
        interactive = false,
        selectedGenes = new Set(),
        onGeneClick = null
    } = options;
    
    // Check if container is in collapsed state
    const isCollapsed = container.classList.contains('heatmap-collapsed');
    const gridClass = isCollapsed ? 'heatmap-grid-collapsed' : 'heatmap-grid';
    
    let html = generateColorbarHTML() + `<div class="${gridClass}">`;
    
    genes.forEach(gene => {
        const expression = expressionData[gene] || 0;
        const uncertainty = uncertaintyData[gene] || 0;
        const expressionClass = getExpressionClass(expression);
        const uncertaintyClass = getUncertaintyClass(uncertainty);
        const isKnockedOut = selectedGenes.has(gene);
        const interactiveClass = interactive ? 'interactive' : '';
        const clickHandler = interactive && onGeneClick ? `onclick="${onGeneClick}('${gene}')"` : '';
        
        html += `
            <div class="heatmap-cell ${expressionClass} ${interactiveClass} ${isKnockedOut ? 'knocked-out' : ''}" 
                 ${clickHandler}
                 data-gene="${gene}"
                 data-expression="${expression.toFixed(2)}"
                 data-uncertainty="${uncertainty.toFixed(2)}"
                 onmouseenter="showHeatmapTooltip(event)"
                 onmouseleave="hideHeatmapTooltip()">
                <div class="heatmap-gene-name">${gene}</div>
                <div class="uncertainty-dot ${uncertaintyClass}"></div>
                ${isKnockedOut ? '<div class="knockout-slash"></div>' : ''}
            </div>
        `;
    });
    
    html += '</div>';
    container.innerHTML = html;
}

/**
 * Show tooltip on heatmap cell hover
 * @param {Event} event - Mouse event
 */
function showHeatmapTooltip(event) {
    const cell = event.currentTarget;
    const gene = cell.dataset.gene;
    const expression = cell.dataset.expression;
    const uncertainty = cell.dataset.uncertainty;
    
    // Remove existing tooltip
    hideHeatmapTooltip();
    
    // Create tooltip
    const tooltip = document.createElement('div');
    tooltip.className = 'heatmap-tooltip show';
    tooltip.innerHTML = `${gene}<br>Mean: ${expression}<br>Std: ${uncertainty}`;
    
    // Position tooltip relative to viewport
    const rect = cell.getBoundingClientRect();
    tooltip.style.position = 'fixed';
    tooltip.style.left = `${rect.left + rect.width / 2}px`;
    tooltip.style.top = `${rect.top - 10}px`;
    tooltip.style.transform = 'translateX(-50%) translateY(-100%)';
    
    document.body.appendChild(tooltip);
}

/**
 * Hide heatmap tooltip
 */
function hideHeatmapTooltip() {
    const existingTooltip = document.querySelector('.heatmap-tooltip');
    if (existingTooltip) {
        existingTooltip.remove();
    }
}

/**
 * Toggle heatmap view between grid and collapsed (single row)
 * @param {string} containerId - ID of the heatmap container
 * @param {Function} rerenderCallback - Function to call to re-render the heatmap
 */
function toggleHeatmapView(containerId, rerenderCallback) {
    const container = document.getElementById(containerId);
    let button;
    
    // Handle different button ID patterns
    if (containerId === 'geneGrid') {
        button = document.getElementById('geneGridCollapseBtn');
    } else if (containerId === 'exploreHeatmap') {
        button = document.getElementById('exploreHeatmapCollapseBtn');
    } else if (containerId.startsWith('exploreHeatmap_')) {
        // Handle numbered explore heatmaps
        const index = containerId.split('_')[1];
        button = document.getElementById(`exploreHeatmapCollapseBtn_${index}`);
    } else {
        button = document.getElementById(containerId.replace('Heatmap', 'CollapseBtn'));
    }
    
    const isCollapsed = container.classList.contains('heatmap-collapsed');
    
    if (isCollapsed) {
        // Expand to grid view
        container.classList.remove('heatmap-collapsed');
        if (button) button.innerHTML = '<i class="fas fa-compress-alt"></i> Collapse';
    } else {
        // Collapse to single row
        container.classList.add('heatmap-collapsed');
        if (button) button.innerHTML = '<i class="fas fa-expand-alt"></i> Expand';
    }
    
    // Re-render the heatmap if callback provided
    if (rerenderCallback) {
        rerenderCallback();
    }
}

/**
 * Show a large loading spinner overlay
 * @param {string} message - Optional loading message
 */
function showSpinner(message = 'Loading...') {
    // Force remove any existing spinner immediately
    hideSpinner();
    
    // Use requestAnimationFrame to ensure DOM is clean before creating new spinner
    requestAnimationFrame(() => {
        const spinner = document.createElement('div');
        spinner.id = 'loadingSpinner';
        spinner.className = 'loading-spinner-overlay';
        spinner.innerHTML = `
            <div class="loading-spinner-content">
                <div class="loading-spinner">
                    <i class="fas fa-spinner fa-spin"></i>
                </div>
                <div class="loading-message">${message}</div>
            </div>
        `;
        
        document.body.appendChild(spinner);
        
        // Force reflow and add show class
        requestAnimationFrame(() => {
            spinner.classList.add('show');
        });
    });
}

/**
 * Hide the loading spinner overlay
 */
function hideSpinner() {
    // Find and remove all spinners (in case there are duplicates)
    const spinners = document.querySelectorAll('#loadingSpinner, .loading-spinner-overlay');
    if (spinners.length > 0) {
        spinners.forEach(spinner => {
            spinner.classList.remove('show');
            // Remove after a short delay to allow fade out
            setTimeout(() => {
                if (spinner.parentNode) {
                    spinner.remove();
                }
            }, 100);
        });
    }
}

// Global instance
const dataService = new DataService();

function toggleSidebar() {
    const sidebar = document.querySelector('.sidebar');
    sidebar.classList.toggle('show');
}

// Close sidebar when clicking outside on mobile
document.addEventListener('click', function(event) {
    const sidebar = document.querySelector('.sidebar');
    const menuBtn = document.querySelector('.mobile-menu-btn');
    
    if (window.innerWidth <= 767.98 && 
        !sidebar.contains(event.target) && 
        !menuBtn.contains(event.target) &&
        sidebar.classList.contains('show')) {
        sidebar.classList.remove('show');
    }
});
