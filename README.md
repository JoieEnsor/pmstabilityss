# pmstabilityss

## Sample Size for Prediction Model Development Based on Prediction Stability

`pmstabilityss` is an R package that calculates the minimum sample size required for developing a clinical prediction model targeting precise individual risk estimates. The package implements the methodology proposed by Riley et al. (2024), which focuses on ensuring that predictions for individual patients have acceptably narrow uncertainty intervals.

## Installation

You can install the development version of pmstabilityss from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("username/pmstabilityss")
```

## Overview

When developing prediction models, it's crucial to have a sufficient sample size not only to prevent overfitting but also to ensure that predictions for individual patients are precise enough to be clinically useful. Traditional sample size calculations often focus on model performance measures like calibration and discrimination, but may not guarantee that individual predictions are sufficiently precise.

This package allows you to:

1. Calculate uncertainty intervals for individual risk predictions at different sample sizes
2. Visualise how prediction stability changes with increasing sample size
3. Determine the sample size needed to achieve target prediction interval widths
4. Compare the uncertainty in classification decisions at different thresholds
5. Analyse prediction stability across different subgroups

## Usage

```r
library(pmstabilityss)

# Basic usage with a dataset containing predictors and a linear predictor
result <- pmstabilityss(data = your_data, 
                        varlist = c("var1", "var2", "var3"), 
                        prevalence = 0.3,
                        lp = "linear_predictor")

# Print results
result$overall_stats

# Plot prediction stability
result$plots$instability_plots

# Sample size needed to achieve specific prediction interval widths
result2 <- pmstabilityss(data = your_data, 
                         varlist = c("var1", "var2", "var3"), 
                         prevalence = 0.3,
                         lp = "linear_predictor",
                         pcutpoints = c(0.3, 0.6, 1),   # Risk categories
                         pciwidth = c(0.1, 0.15, 0.2))  # Target widths for each category

# Classification instability at threshold 0.5
result3 <- pmstabilityss(data = your_data, 
                         varlist = c("var1", "var2", "var3"), 
                         prevalence = 0.3,
                         lp = "linear_predictor",
                         threshold = 0.5)

# Plot classification instability
result3$plots$classification_plots
```

## Key Parameters

- **data**: Your dataset containing predictor variables
- **varlist**: Character vector of predictor variable names
- **prevalence**: Expected prevalence of the outcome
- **lp**: Name of the linear predictor variable (optional)
- **pcutpoints** and **pciwidth**: Specify target prediction interval widths for different probability ranges
- **threshold**: Risk threshold for evaluating classification stability
- **pmss**: Set to TRUE to use the pmsampsize package for comparison
- **n**: Additional sample sizes to evaluate

## Output

The function returns an object of class "pmstabilityss" containing:

- **overall_stats**: Summary statistics of prediction interval widths for each sample size
- **threshold_stats**: Summary of misclassification probabilities (if threshold specified)
- **cuts_stats**: Summary by probability categories (if pciwidth specified)
- **plots**: Visualisations of prediction stability and classification uncertainty

## Citation

Riley, R. D., Collins, G. S., Whittle, R., Archer, L., Snell, K. I., Dhiman, P., ... & Ensor, J. (2024). Sample size for developing a prediction model with a binary outcome: targeting precise individual risk estimates to improve clinical decisions and fairness. *arXiv preprint arXiv:2407.09293*.

## Authors

Joie Ensor, University of Birmingham (j.ensor@bham.ac.uk)
