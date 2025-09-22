// Global variables
let datasets = [];
let genes = [];
let selectedTarget = null;
let knockoutBudget = 2;

// Load data and initialize page
document.addEventListener('DOMContentLoaded', function() {
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
    }
}

// Load genes using DataService
async function loadGenes() {
    try {
        genes = await dataService.loadGenes();
        genes.sort(); // Sort alphabetically
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
        const tissueOrgans = await dataService.getUniqueValues('tissue_or_organ');

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
}

// Update UMAP visualization based on filters
function updateOptimizeVisualization() {
    const filters = {
        tcellSubtype: document.getElementById('optimizeTcellSubtype').value,
        donorType: document.getElementById('optimizeDonorType').value,
        timePoint: document.getElementById('optimizeTimePoint').value,
        tissueOrgan: document.getElementById('optimizeTissueOrgan').value
    };

    // Get color by option
    const colorBy = document.getElementById('optimizeColorBy').value;

    // Generate image filename based on filters
    const imageName = generateOptimizeImageName(filters, colorBy);
    // Add timestamp to force new image generation
    const timestamp = Date.now();
    const imagePath = `images/optimize/umaps/${imageName}?t=${timestamp}`;
    
    document.getElementById('optimizeUmapPlot').src = imagePath;
    
    // Reset target selection when filters change
    resetTargetSelection();
}

// Generate optimize image filename based on filters
function generateOptimizeImageName(filters, colorBy = 'tcell_subtype') {
    const parts = [];
    
    if (filters.tcellSubtype !== 'all') parts.push(`tcell-${filters.tcellSubtype.replace(/\s+/g, '-').toLowerCase()}`);
    if (filters.donorType !== 'all') parts.push(`donor-${filters.donorType.replace(/\s+/g, '-').toLowerCase()}`);
    if (filters.timePoint !== 'all') parts.push(`time-${filters.timePoint.replace(/\s+/g, '-').toLowerCase()}`);
    if (filters.tissueOrgan !== 'all') parts.push(`tissue-${filters.tissueOrgan.replace(/\s+/g, '-').toLowerCase()}`);
    
    // Add color by parameter
    parts.push(`color-${colorBy.replace(/_/g, '-')}`);
    
    return parts.length > 0 ? `${parts.join('_')}.png` : 'default.png';
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
    document.getElementById('budgetDisplay2').textContent = knockoutBudget;
}

// Handle target area selection on UMAP
function selectTargetArea(event) {
    const img = event.target;
    const rect = img.getBoundingClientRect();
    
    // Calculate relative coordinates (0-1 range)
    const x = (event.clientX - rect.left) / rect.width;
    const y = (event.clientY - rect.top) / rect.height;
    
    // Convert to UMAP coordinates (approximate range -10 to 10)
    const umapX = (x - 0.5) * 20;
    const umapY = (0.5 - y) * 20; // Flip Y axis
    
    selectedTarget = { x: umapX, y: umapY, pixelX: x, pixelY: y };
    
    // Update UI
    showTargetOverlay(x * rect.width, y * rect.height);
    updateTargetInfo(umapX, umapY);
}

// Show target overlay on the UMAP
function showTargetOverlay(pixelX, pixelY) {
    const overlay = document.getElementById('targetOverlay');
    const img = document.getElementById('optimizeUmapPlot');
    const rect = img.getBoundingClientRect();
    
    overlay.style.display = 'block';
    overlay.style.left = `${pixelX - 15}px`; // Center the 30px circle
    overlay.style.top = `${pixelY - 15}px`;
}

// Update target info display
function updateTargetInfo(umapX, umapY) {
    document.getElementById('targetCoords').textContent = `(${umapX.toFixed(1)}, ${umapY.toFixed(1)})`;
    document.getElementById('targetInfo').style.display = 'block';
}


// Reset target selection
function resetTargetSelection() {
    selectedTarget = null;
    document.getElementById('targetOverlay').style.display = 'none';
    document.getElementById('targetInfo').style.display = 'none';
}

// Run optimization algorithm
function runOptimization() {
    if (!selectedTarget) {
        alert('Please select a target area on the UMAP first.');
        return;
    }
    
    // Show loading state
    const button = document.getElementById('optimizeButton');
    const originalText = button.innerHTML;
    button.innerHTML = '<i class="fas fa-spinner fa-spin"></i> Optimizing...';
    button.disabled = true;
    
    // Simulate optimization process
    setTimeout(() => {
        const results = generateOptimizationResults();
        displayResults(results);
        
        // Reset button
        button.innerHTML = originalText;
        button.disabled = false;
        
        // Scroll to results
        document.getElementById('resultsCard').scrollIntoView({ 
            behavior: 'smooth', 
            block: 'start' 
        });
    }, 2000); // 2 second delay to simulate computation
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
    
    // Update results UMAP image
    updateResultsUmap();
    
    // Show result overlays
    showResultOverlays(results, resultColors);
    
    document.getElementById('resultsCard').style.display = 'block';
}

// Update the UMAP image in results to match current filters
function updateResultsUmap() {
    const currentUmapSrc = document.getElementById('optimizeUmapPlot').src;
    const resultsImg = document.getElementById('resultsUmapPlot');
    resultsImg.src = currentUmapSrc;
    
    // Wait for image to load before positioning overlays
    resultsImg.onload = function() {
        // Small delay to ensure DOM is settled
        setTimeout(() => {
            if (selectedTarget) {
                const rect = resultsImg.getBoundingClientRect();
                const overlay = document.getElementById('originalTargetOverlay');
                
                overlay.style.display = 'block';
                overlay.style.left = `${selectedTarget.pixelX * rect.width - 15}px`;
                overlay.style.top = `${selectedTarget.pixelY * rect.height - 15}px`;
            }
        }, 50);
    };
}

// Show result overlays on the UMAP
function showResultOverlays(results, colors) {
    const img = document.getElementById('resultsUmapPlot');
    
    // Ensure image is loaded before positioning overlays
    const positionOverlays = () => {
        const rect = img.getBoundingClientRect();
        
        // Only proceed if we have valid dimensions
        if (rect.width === 0 || rect.height === 0) {
            setTimeout(positionOverlays, 100);
            return;
        }
        
        results.forEach((result, index) => {
            // Generate random positions within the plot boundaries (0.1 to 0.9 range to avoid edges)
            const resultX = 0.1 + Math.random() * 0.8;
            const resultY = 0.1 + Math.random() * 0.8;
            
            const overlay = document.getElementById(`result${index + 1}Overlay`);
            
            overlay.style.display = 'block';
            overlay.style.left = `${resultX * rect.width - 12}px`;
            overlay.style.top = `${resultY * rect.height - 12}px`;
            overlay.style.backgroundColor = colors[index];
            overlay.style.borderColor = colors[index];
        });
    };
    
    // Wait a bit to ensure the image and layout are ready
    setTimeout(positionOverlays, 100);
}
