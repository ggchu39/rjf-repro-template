# Workflow: Gamma nested bootstrap CV template

## What this template is for
This repository provides a reusable analysis pattern for positive, right-skewed outcomes modelled using Gamma regression with a log link.

It currently focuses on:
- Outer loop: bootstrap resampling with LOOCV-style folds to estimate prediction error distributions.
- Inner loop: candidate-model comparison using AICc (information-theoretic selection) within each fold.

## Why nested selection
Selecting a model using all data and then evaluating it on the same data can overstate performance.
Nested selection keeps model choice separated from performance estimation.

## Why fold-wise scaling
If predictors are standardized using the full dataset, information from the held-out fold can leak into training.
This template standardizes numeric predictors using the outer training fold only, then applies the same transform to the outer test fold.

## Optional diagnostics (may be added/extended)
- Convergence monitoring (inner fits and outer refits)
- Error summaries (MAE/RMSE distributions, tail behaviour, worst-case checks)

## Trade-offs / limitations
- Computationally expensive (B bootstraps × n folds × M models).
- AICc focuses on relative fit and parsimony; it is not a direct measure of predictive accuracy.
- If many fits fail to converge, consider simplifying candidate models or adjusting optimizer controls.
