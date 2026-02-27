# DECISIONS log

## 2025-11-15
- Adopted a Gamma(log) modelling pipeline for positive, right-skewed latency-style outcomes.
- Replaced simple LOOCV with a nested bootstrap CV design (outer bootstrap + LOOCV-style folds) to better separate model selection from performance estimation.
- Implemented inner-loop AICc model selection to prioritise parsimonious candidate models under small-sample conditions.
- Added fold-wise scaling for numeric predictors to reduce information leakage across folds.
- Added convergence monitoring (inner fits and outer refits) and guardrails for invalid predictions (default: return NA; optional clipping).

## 2026-02-27
- Published a public-facing template repository focusing on the Gamma(log) nested bootstrap CV “engine” and documentation; project-specific data, model lists, and outputs remain private.
