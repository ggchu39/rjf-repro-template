# Workflow: Negative Binomial paired nested LOO-in-bootstrap CV (template)

## What this template is for
This template evaluates **count outcomes** using **Negative Binomial regression** (MASS::glm.nb)
and estimates out-of-sample performance using a paired resampling design.

It is designed to be **data-agnostic**: you provide a data.frame and a named list of formulas.

## Core idea (leakage-proof evaluation)
For each bootstrap iteration:
1. Draw a bootstrap sample (same resample used for every candidate model → “paired”).
2. For each model separately, run **leave-one-out (LOO)** inside that bootstrap sample.
3. Compute MAE/RMSE from the out-of-fold predictions.

Because each prediction is made on a held-out case, performance is evaluated out-of-fold
rather than on the same data used to fit that fold.

## Outputs (source of truth)
`cvb_nb_single_nestedLOO_paired()` writes one RDS:

- `SingleCVB_Paired_<region_label>.rds`

This RDS is the **source of truth** and contains:
- `mae` and `rmse`: B × M matrices (bootstraps × models)
- `summary_raw`: untrimmed summaries
- `convergence`: fold-level success counts and success rate
- `meta`: B, seed, outcome_var, region_label

All tables are derived from this saved RDS to ensure traceability.

## Reporting vs computation
The engine saves untrimmed results.
For reporting only, `robust_summarize()` trims catastrophic values (non-finite or > cutoffs)
and reports how many were dropped. The underlying RDS is not modified.

## Manuscript-style tables
- `make_table_312()`: per-model MAE/RMSE with bootstrap percentile CIs + success rate
- `make_table_313()`: paired model comparisons (Δ, 95% CI, win%, N_pairs)
  - Deltas are computed **on paired rows only** (same bootstrap iterations where both models are finite).

## Practical notes
- Runtime scales as: B × n × M (bootstraps × LOO folds × models).
- Use small B for demos; use manuscript B for final runs.
- If fold-level success rates are low, simplify model formulas or reduce separation/instability.
