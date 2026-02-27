# DECISIONS (index)

This file is a short index of key implementation decisions for the CV templates in this repository.

## Per-template decision logs
- Gamma: [docs/Gamma_DECISIONS.md](docs/Gamma_DECISIONS.md)
- Logistic: [docs/Logit_DECISIONS.md](docs/Logit_DECISIONS.md)
- Negative Binomial: [docs/NB_DECISIONS.md](docs/NB_DECISIONS.md)

## Global decisions (applies across templates)
- **Leakage control:** performance is always estimated out-of-sample (no training-on-test reuse).
- **Reproducibility:** explicit random seeds are passed and recorded; outputs are deterministic given the same inputs.
- **Source of truth:** model-evaluation results are saved once (RDS) and all tables are derived from that saved object.
- **Paired resampling (NB):** all candidate models are evaluated on the *same* bootstrap resamples to enable paired comparisons.
- **Reporting-only trimming:** any “catastrophic value” trimming is applied only for reporting summaries; the original stored results remain unchanged.
