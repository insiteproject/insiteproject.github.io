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
        
        // Count unique records (repositories) instead of data objects
        const uniqueRecords = new Set(datasets.map(d => d.repo_id)).size;
        
        return {
            recordCount: uniqueRecords,
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
        const filteredDatasets = await this.filterDatasets(filters);
        const genes = await this.loadGenes();
        
        if (filteredDatasets.length === 0) {
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
        filteredDatasets.forEach((dataset) => {
            const weight = dataset.n_obs || 1; // Use cell count as weight
            totalWeight += weight;
            
            if (dataset.mean_expression_profile && dataset.mean_expression_profile.length === genes.length) {
                genes.forEach((gene, index) => {
                    const profile = dataset.mean_expression_profile[index];
                    if (profile && Array.isArray(profile) && profile.length === 2) {
                        weightedMeans[gene] += profile[0] * weight; // mean is first element
                        weightedStds[gene] += profile[1] * weight;  // std is second element
                    }
                });
            }
        });
        
        // Normalize by total weight
        const result = {};
        genes.forEach(gene => {
            result[gene] = {
                mean: totalWeight > 0 ? weightedMeans[gene] / totalWeight : 0,
                std: totalWeight > 0 ? weightedStds[gene] / totalWeight : 0
            };
        });
        
        return result;
    }
}

/**
 * Generate colorbar HTML for heatmaps
 * @returns {string} HTML string for colorbar
 */
function generateColorbarHTML() {
    return `
        <div class="heatmap-colorbars">
            <div class="colorbar-container">
                <div class="colorbar-label">Expression (Mean)</div>
                <div class="colorbar expression-colorbar">
                    <div class="colorbar-item expr-low">0-2</div>
                    <div class="colorbar-item expr-low-med">2-4</div>
                    <div class="colorbar-item expr-medium">4-6</div>
                    <div class="colorbar-item expr-med-high">6-8</div>
                    <div class="colorbar-item expr-high">8-10</div>
                </div>
            </div>
            <div class="colorbar-container">
                <div class="colorbar-label">Uncertainty (Std)</div>
                <div class="colorbar uncertainty-colorbar">
                    <div class="colorbar-item uncert-very-low">0-0.4</div>
                    <div class="colorbar-item uncert-low">0.4-0.8</div>
                    <div class="colorbar-item uncert-medium">0.8-1.2</div>
                    <div class="colorbar-item uncert-high">1.2-1.6</div>
                    <div class="colorbar-item uncert-very-high">1.6-2.0</div>
                </div>
            </div>
        </div>`;
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
