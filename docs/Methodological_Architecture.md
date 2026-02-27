
1. Outcomes & model families

- Count Outcomes: Negative Binomial GLMs.
- Skewed Latency: Gamma GLMs (log link).
- Binary Success: Logistic GLMs.

2. Candidate Modelling & Collinearity Control
A Priori Specification: Candidate models defined as biologically plausible combinations (2–6 predictors), including additive and selected two-way interaction structures.
Collinearity Screening: Rigorous assessment using Variance Inflation Factors (VIF) and Adjusted Generalized VIF (GVIF) for multi-parameter terms.
Heuristic Thresholds: Strict monitoring of GVIF values (threshold ~2–3) alongside model stability and diagnostic performance to ensure robust parameter estimation.

3. Information-Theoretic Model Ranking
Selection Framework: Models ranked using AICc (Akaike Information Criterion corrected for small samples).
Shortlisting: Retention of an information-theoretic shortlist (ΔAICc ≤ 4) for further evaluation.
Prioritisation: Final inference prioritises parsimony, biological plausibility, and diagnostic integrity.

4. Diagnostics & fit checks

- DHARMa simulated-residual diagnostics for count/Gamma models (dispersion, zero inflation, residual patterns).
- Logistic models assessed with discrimination (AUC) and calibration (Hosmer–Lemeshow + graphical calibration).
- Standard diagnostic also applied.

5. Leakage-proof validation (core quality signal)

- Nested resampling to separate model choice from performance estimation (avoid information leakage).
- Two-stage validation:
  (i)  Multi-model nested resampling for model selection + out-of-fold performance.
  (ii) Single-model resampling for uncertainty / confidence intervals of performance metrics.

6. Paired resampling for fair comparisons
When comparing multiple fixed models, the workflow uses paired bootstrap resamples across models so head-to-head differences are computed on the same resamples.

7. Robust reporting

- “Catastrophic” or non-finite resampling errors are trimmed for reporting only, while preserving the source-of-truth saved object.
- All reporting tables are derived from the saved outputs (no refitting during reporting).

8. Sensitivity analysis
Turnover robustness tested by swapping DOPAC/DA ↔ HVA/DA in otherwise identical model backbones and re-running the same validation pipeline.

9. Reproducibility
Deterministic seeds, consistent naming, and a “source-of-truth RDS → derived tables” structure for traceability.

