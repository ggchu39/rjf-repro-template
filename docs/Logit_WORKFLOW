# Workflow: Logistic nested bootstrap CV template

## What this template is for
This repository provides a reusable analysis pattern for a **binary outcome** (0/1)
modelled using **logistic regression**.

It implements:
- **Outer loop:** bootstrap resampling to estimate variability in out-of-sample performance.
- **Inner loop (within each bootstrap):** leave-one-out folds where candidate models are compared
  using **AICc**, and the best model is used to predict the held-out case.

## Why nested selection
Selecting a model on the full dataset and evaluating it on the same dataset can
overstate performance.
Nested selection keeps **model choice** (inner loop) separated from **performance estimation**
(outer loop).

## Metrics returned
Across bootstraps:
- **AUC** (discrimination; computed only when both classes appear in out-of-fold pairs)
- **Brier score** (probability accuracy)
- **Log loss** (probabilistic scoring rule)
- **Calibration intercept & slope** (from recalibration on out-of-fold predictions)
- **Hosmer–Lemeshow p-value** (with adaptive bin count, when feasible)

Also returned:
- **Model selection frequency** (most commonly selected model within each bootstrap)
- **Convergence matrix** (inner-loop fit success per model per bootstrap)
- **Per-bootstrap out-of-fold predictions** (for optional plots)

## Trade-offs / limitations
- Computationally expensive (B bootstraps × n folds × M models).
- AICc is a relative-fit criterion, not a direct measure of predictive accuracy.
- With small samples, separation or near-separation can cause non-convergence; simplify
  candidate models if convergence rates are low.
- Hosmer–Lemeshow is sensitive and may be unstable in small samples; treat it as a rough diagnostic.
