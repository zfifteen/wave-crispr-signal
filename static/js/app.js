// JavaScript functions for Z Framework web application

// Global variables
let loadingSpinner = null;

// Initialize application
document.addEventListener('DOMContentLoaded', function() {
    initializeTooltips();
    initializeLoadingSpinner();
    setupFormValidation();
    animateElements();
});

// Initialize Bootstrap tooltips
function initializeTooltips() {
    const tooltipTriggerList = [].slice.call(document.querySelectorAll('[data-bs-toggle="tooltip"]'));
    tooltipTriggerList.map(function (tooltipTriggerEl) {
        return new bootstrap.Tooltip(tooltipTriggerEl);
    });
}

// Initialize loading spinner
function initializeLoadingSpinner() {
    loadingSpinner = document.createElement('div');
    loadingSpinner.className = 'loading-spinner';
    loadingSpinner.innerHTML = `
        <div class="spinner-border text-primary" role="status">
            <span class="visually-hidden">Loading...</span>
        </div>
        <div class="mt-2 text-center">
            <small class="text-muted">Calculating Z Framework metrics...</small>
        </div>
    `;
    document.body.appendChild(loadingSpinner);
}

// Setup form validation
function setupFormValidation() {
    const forms = document.querySelectorAll('form[method="POST"]');
    forms.forEach(form => {
        form.addEventListener('submit', function(event) {
            if (!validateDNASequence(form)) {
                event.preventDefault();
                event.stopPropagation();
            } else {
                showLoadingSpinner();
            }
        });
    });
}

// Validate DNA sequence
function validateDNASequence(form) {
    const sequenceField = form.querySelector('[name="sequence"]');
    if (!sequenceField) return true;
    
    const sequence = sequenceField.value.trim().toUpperCase();
    const validBases = /^[ATCG]+$/;
    
    // Clear previous validation
    sequenceField.classList.remove('is-invalid');
    const existingFeedback = sequenceField.parentNode.querySelector('.invalid-feedback');
    if (existingFeedback) {
        existingFeedback.remove();
    }
    
    // Validate sequence
    if (sequence.length < 10) {
        showValidationError(sequenceField, 'Sequence must be at least 10 bases long.');
        return false;
    }
    
    if (sequence.length > 1000) {
        showValidationError(sequenceField, 'Sequence must be no more than 1000 bases long.');
        return false;
    }
    
    if (!validBases.test(sequence)) {
        showValidationError(sequenceField, 'Sequence must contain only A, T, C, G bases.');
        return false;
    }
    
    return true;
}

// Show validation error
function showValidationError(field, message) {
    field.classList.add('is-invalid');
    const feedback = document.createElement('div');
    feedback.className = 'invalid-feedback';
    feedback.textContent = message;
    field.parentNode.appendChild(feedback);
}

// Show loading spinner
function showLoadingSpinner() {
    if (loadingSpinner) {
        loadingSpinner.style.display = 'block';
    }
}

// Hide loading spinner
function hideLoadingSpinner() {
    if (loadingSpinner) {
        loadingSpinner.style.display = 'none';
    }
}

// Animate elements on page load
function animateElements() {
    const cards = document.querySelectorAll('.card');
    cards.forEach((card, index) => {
        setTimeout(() => {
            card.classList.add('fade-in');
        }, index * 100);
    });
}

// Format DNA sequence with base coloring
function formatDNASequence(sequence, elementId) {
    const element = document.getElementById(elementId);
    if (!element) return;
    
    let formattedSequence = '';
    for (let base of sequence) {
        const className = `base-${base.toLowerCase()}`;
        formattedSequence += `<span class="${className}">${base}</span>`;
    }
    
    element.innerHTML = formattedSequence;
}

// Copy text to clipboard
function copyToClipboard(text, showToast = true) {
    navigator.clipboard.writeText(text).then(function() {
        if (showToast) {
            showToast('Text copied to clipboard!', 'success');
        }
    }).catch(function() {
        // Fallback for older browsers
        const textArea = document.createElement('textarea');
        textArea.value = text;
        document.body.appendChild(textArea);
        textArea.select();
        document.execCommand('copy');
        document.body.removeChild(textArea);
        
        if (showToast) {
            showToast('Text copied to clipboard!', 'success');
        }
    });
}

// Show toast notification
function showToast(message, type = 'info') {
    const toast = document.createElement('div');
    toast.className = `toast position-fixed top-0 end-0 m-3`;
    toast.setAttribute('role', 'alert');
    toast.innerHTML = `
        <div class="toast-header">
            <i class="bi bi-${type === 'success' ? 'check-circle' : 'info-circle'} text-${type} me-2"></i>
            <strong class="me-auto">${type === 'success' ? 'Success' : 'Info'}</strong>
            <button type="button" class="btn-close" data-bs-dismiss="toast"></button>
        </div>
        <div class="toast-body">
            ${message}
        </div>
    `;
    
    document.body.appendChild(toast);
    const bsToast = new bootstrap.Toast(toast);
    bsToast.show();
    
    // Remove toast element after it's hidden
    toast.addEventListener('hidden.bs.toast', () => {
        document.body.removeChild(toast);
    });
}

// Download data as JSON file
function downloadJSON(data, filename) {
    const dataStr = "data:text/json;charset=utf-8," + encodeURIComponent(JSON.stringify(data, null, 2));
    const downloadAnchorNode = document.createElement('a');
    downloadAnchorNode.setAttribute("href", dataStr);
    downloadAnchorNode.setAttribute("download", filename);
    document.body.appendChild(downloadAnchorNode);
    downloadAnchorNode.click();
    downloadAnchorNode.remove();
}

// Calculate base composition
function calculateBaseComposition(sequence) {
    const composition = { A: 0, T: 0, C: 0, G: 0 };
    for (let base of sequence.toUpperCase()) {
        if (composition.hasOwnProperty(base)) {
            composition[base]++;
        }
    }
    
    const total = sequence.length;
    const percentages = {};
    for (let base in composition) {
        percentages[base] = (composition[base] / total * 100).toFixed(1);
    }
    
    return {
        counts: composition,
        percentages: percentages,
        gcContent: ((composition.G + composition.C) / total * 100).toFixed(1)
    };
}

// Format number with specific precision
function formatNumber(number, precision = 6) {
    if (typeof number === 'number') {
        return number.toFixed(precision);
    }
    return number;
}

// Generate random DNA sequence
function generateRandomDNA(length, gcContent = 0.5) {
    const gcBases = ['G', 'C'];
    const atBases = ['A', 'T'];
    let sequence = '';
    
    for (let i = 0; i < length; i++) {
        if (Math.random() < gcContent) {
            sequence += gcBases[Math.floor(Math.random() * gcBases.length)];
        } else {
            sequence += atBases[Math.floor(Math.random() * atBases.length)];
        }
    }
    
    return sequence;
}

// Smooth scroll to element
function scrollToElement(elementId) {
    const element = document.getElementById(elementId);
    if (element) {
        element.scrollIntoView({
            behavior: 'smooth',
            block: 'start'
        });
    }
}

// API call wrapper
async function callAnalysisAPI(sequence, parameters = {}) {
    const defaultParams = {
        precision: 50,
        window_size: 10,
        num_perturbations: 100,
        perturbation_rate: 0.1
    };
    
    const requestData = {
        sequence: sequence,
        ...defaultParams,
        ...parameters
    };
    
    try {
        const response = await fetch('/api/analyze', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify(requestData)
        });
        
        if (!response.ok) {
            throw new Error(`HTTP error! status: ${response.status}`);
        }
        
        return await response.json();
    } catch (error) {
        console.error('API call failed:', error);
        throw error;
    }
}

// Chart utilities
const ChartUtils = {
    // Common chart options
    getCommonOptions: function(title) {
        return {
            responsive: true,
            maintainAspectRatio: false,
            plugins: {
                title: {
                    display: true,
                    text: title,
                    font: {
                        size: 16,
                        weight: 'bold'
                    }
                },
                legend: {
                    position: 'top'
                }
            },
            scales: {
                y: {
                    beginAtZero: false,
                    grid: {
                        color: 'rgba(0, 0, 0, 0.1)'
                    }
                },
                x: {
                    grid: {
                        color: 'rgba(0, 0, 0, 0.1)'
                    }
                }
            }
        };
    },
    
    // Color palette
    colors: {
        primary: 'rgb(75, 192, 192)',
        secondary: 'rgb(54, 162, 235)',
        success: 'rgb(75, 192, 75)',
        warning: 'rgb(255, 206, 86)',
        danger: 'rgb(255, 99, 132)',
        info: 'rgb(54, 162, 235)'
    }
};

// Hide loading spinner when page loads
window.addEventListener('load', hideLoadingSpinner);