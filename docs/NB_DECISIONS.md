## NB CV template (public)

- Chose a paired bootstrap design (same resamples for all candidate models) to support fair per-bootstrap comparisons.
- Used LOO within each bootstrap sample to generate out-of-fold predictions for MAE/RMSE.
- Saved a single RDS per run as the “source of truth”; all reporting tables are derived from that artifact.
- Kept robust trimming as reporting-only (RDS remains unmodified) to preserve traceability.
- Public repo contains generic engines + toy demos; project-specific model lists and real-data scripts are kept private.
