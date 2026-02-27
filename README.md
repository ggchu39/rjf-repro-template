**A Robust Framework for Multi-Model Inference and Leakage-Proof Nested Validation (Behavioural Neuroscience)**
This repository provides portfolio-grade R templates for reproducible modelling and leakage-proof nested bootstrap cross-validation, tailored to common behavioural neuroscience outcomes:
- Counts → Negative Binomial GLMs
- Skewed latencies → Gamma GLMs (log link)
- Binary success → Logistic GLMs

It includes reporting workflows (source-of-truth saved objects → derived tables) and paired resampling utilities for fair head-to-head model comparisons.

**Note:** This is a methodological template illustrated with simulated/toy data. It intentionally excludes raw project data, manuscript-specific model specifications, and results to respect privacy and publication embargoes while documenting reusable methodology and authorship.

To use the template, provide your own dataframe with:
- a positive outcome variable (e.g., latency-like measure)
- predictors referenced by the candidate model formulas
