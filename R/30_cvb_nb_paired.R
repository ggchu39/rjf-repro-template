# R/cvb_nb_paired.R
# Negative Binomial (count) paired nested LOO-in-bootstrap CV
# Public template: generic + data-agnostic
# - Outer loop: bootstrap resampling
# - Inner loop: leave-one-out within each bootstrap sample (per model, independently)
# - Outputs: RDS (source of truth), robust summaries, Table 312/313 builders

#' Safe Negative Binomial fit with warning capture
#'
#' Fits a negative binomial GLM via MASS::glm.nb and captures warnings.
#' The returned object is designed to be checked by `nb_ok()`.
#'
#' @param form model formula.
#' @param df data.frame containing outcome + predictors referenced in `form`.
#' @param control glm.control passed to glm.nb.
#' @param maxit_nb maximum iterations for NB fitting (glm.nb `maxit`).
#' @return A list with elements:
#'   \itemize{
#'     \item fit: fitted model object or "try-error"
#'     \item warns: character vector of warning messages captured during fitting
#'   }
nb_fit_safely <- function(form, df,
                          control = stats::glm.control(maxit = 200),
                          maxit_nb = 300) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for glm.nb().")
  }

  warns <- character(0)
  fit <- withCallingHandlers(
    try(MASS::glm.nb(form, data = df, control = control, maxit = maxit_nb,
                     na.action = stats::na.exclude),
        silent = TRUE),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(fit = fit, warns = warns)
}

#' Check whether an NB fit is usable
#'
#' Flags failed fits, non-convergence, and common instability warnings.
#'
#' @param fit_obj list returned by `nb_fit_safely()`.
#' @return Logical scalar. TRUE if fit is acceptable for prediction; otherwise FALSE.
nb_ok <- function(fit_obj) {
  if (!is.list(fit_obj) || is.null(fit_obj$fit)) return(FALSE)
  if (inherits(fit_obj$fit, "try-error")) return(FALSE)

  # Common instability warnings (adjust to taste; keep conservative for public template)
  if (length(fit_obj$warns)) {
    bad <- grepl("alternation limit reached|iteration limit reached|NaNs produced",
                 fit_obj$warns, ignore.case = TRUE)
    if (any(bad)) return(FALSE)
  }

  if (!is.null(fit_obj$fit$converged) && isFALSE(fit_obj$fit$converged)) return(FALSE)

  co <- try(stats::coef(fit_obj$fit), silent = TRUE)
  if (inherits(co, "try-error") || any(!is.finite(co))) return(FALSE)

  TRUE
}

#' Leave-one-out prediction within a single bootstrap sample for one NB model
#'
#' For a given bootstrap sample and a single model formula, performs LOO fits and
#' produces out-of-fold predictions, then computes MAE and RMSE.
#'
#' @param df_boot data.frame bootstrap sample.
#' @param form a single model formula.
#' @param outcome_var character name of the outcome column (count).
#' @param control glm.control for fitting.
#' @param maxit_nb maximum iterations for NB fitting.
#' @return List with elements:
#'   \itemize{
#'     \item mae: MAE across valid folds
#'     \item rmse: RMSE across valid folds
#'     \item ok_fits: number of successful fold fits (valid prediction produced)
#'     \item total_fits: number of LOO folds attempted (= nrow(df_boot))
#'   }
nested_loo_one_model <- function(df_boot, form, outcome_var,
                                 control = stats::glm.control(maxit = 200),
                                 maxit_nb = 300) {
  if (!requireNamespace("Metrics", quietly = TRUE)) {
    stop("Package 'Metrics' is required for MAE/RMSE.")
  }

  n_b <- nrow(df_boot)
  y   <- df_boot[[outcome_var]]
  pred <- rep(NA_real_, n_b)
  ok_count <- 0L

  for (i in seq_len(n_b)) {
    tr <- df_boot[-i, , drop = FALSE]
    te <- df_boot[ i, , drop = FALSE]

    fit <- nb_fit_safely(form, tr, control = control, maxit_nb = maxit_nb)
    if (nb_ok(fit)) {
      pr <- try(stats::predict(fit$fit, newdata = te, type = "response"),
                silent = TRUE)
      if (!inherits(pr, "try-error") && length(pr) == 1L && is.finite(pr)) {
        pred[i] <- as.numeric(pr)
        ok_count <- ok_count + 1L
      }
    }
  }

  good <- is.finite(pred) & is.finite(y)
  list(
    mae  = if (!any(good)) NA_real_ else Metrics::mae(y[good],  pred[good]),
    rmse = if (!any(good)) NA_real_ else Metrics::rmse(y[good], pred[good]),
    ok_fits = ok_count,
    total_fits = n_b
  )
}

#' Paired nested LOO-in-bootstrap CV for a fixed NB model set
#'
#' Evaluates each candidate model independently using the SAME bootstrap resamples
#' ("paired" across models). The saved RDS contains the full MAE/RMSE matrices and
#' is treated as the source of truth for all downstream tables.
#'
#' @param data data.frame containing outcome + predictors.
#' @param model_formulas named list of formulas. If unnamed, names are created.
#' @param region_label character tag used for output file prefixes.
#' @param outcome_var character name of the outcome column.
#' @param B integer number of bootstrap iterations.
#' @param seed integer RNG seed.
#' @param out_dir output directory (default ".").
#' @return Invisibly returns a list with rds_path and basic in-memory summary.
cvb_nb_single_nestedLOO_paired <- function(data, model_formulas, region_label,
                                          outcome_var,
                                          B = 2000, #example from my manuscript
                                          seed = 190825,#example from my manuscript
                                          out_dir = ".") {
  stopifnot(is.data.frame(data))
  stopifnot(length(model_formulas) >= 1L)
  stopifnot(is.character(region_label), length(region_label) == 1L)
  stopifnot(is.character(outcome_var), length(outcome_var) == 1L)
  stopifnot(is.numeric(B), length(B) == 1L, B >= 1L)
  stopifnot(is.numeric(seed), length(seed) == 1L)

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  labs <- names(model_formulas)
  if (is.null(labs) || any(labs == "")) {
    names(model_formulas) <- paste0("Model ", seq_along(model_formulas))
    labs <- names(model_formulas)
  }

  set.seed(seed)
  n <- nrow(data)
  M <- length(model_formulas)

  mae_mat  <- matrix(NA_real_, nrow = B, ncol = M, dimnames = list(NULL, labs))
  rmse_mat <- matrix(NA_real_, nrow = B, ncol = M, dimnames = list(NULL, labs))
  ok_vec   <- integer(M)
  tot_vec  <- integer(M)

  for (b in seq_len(B)) {
    idx  <- sample.int(n, replace = TRUE)
    boot <- data[idx, , drop = FALSE]

    for (m in seq_len(M)) {
      res <- nested_loo_one_model(boot, model_formulas[[m]], outcome_var = outcome_var)
      mae_mat[b, m]  <- res$mae
      rmse_mat[b, m] <- res$rmse
      ok_vec[m]  <- ok_vec[m]  + res$ok_fits
      tot_vec[m] <- tot_vec[m] + res$total_fits
    }
  }

  # Raw summaries (no trimming here; trimming is reporting-only)
  sum_tbl <- data.frame(
    Model     = labs,
    Mean_MAE  = colMeans(mae_mat,  na.rm = TRUE),
    Mean_RMSE = colMeans(rmse_mat, na.rm = TRUE),
    Max_MAE   = apply(mae_mat,  2, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)),
    Max_RMSE  = apply(rmse_mat, 2, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
  )

  conv_tbl <- data.frame(
    Model = labs,
    SuccessfulFits = ok_vec,
    TotalFitsAttempted = tot_vec,
    SuccessRate_FoldLevel = ok_vec / pmax(1L, tot_vec)
  )

  rds_path <- file.path(out_dir, paste0("SingleCVB_Paired_", region_label, ".rds"))
  saveRDS(list(
    mae = mae_mat,
    rmse = rmse_mat,
    summary_raw = sum_tbl,
    convergence = conv_tbl,
    model_order = labs,
    meta = list(B = B, seed = seed, outcome_var = outcome_var, region_label = region_label)
  ), file = rds_path)

  invisible(list(rds_path = rds_path, summary_raw = sum_tbl))
}

#' Robust reporting summary (trim catastrophic values only for reporting)
#'
#' Applies a cutoff to MAE/RMSE matrices stored in the RDS and writes a *_ROBUST.csv.
#' This does NOT change the source-of-truth RDS.
#'
#' @param r either the list read from the RDS, or a path to the RDS file.
#' @param region_label label used in output filename.
#' @param cut_mae values > cut_mae (or non-finite) are set to NA for reporting.
#' @param cut_rmse values > cut_rmse (or non-finite) are set to NA for reporting.
#' @param out_dir output directory.
#' @return List with trimmed matrices and the robust summary table.
robust_summarize <- function(r, region_label,
                             cut_mae = 1e3, cut_rmse = 1e3,
                             out_dir = ".") {
  if (is.character(r)) r <- readRDS(r)
  stopifnot(is.list(r), !is.null(r$mae), !is.null(r$rmse))

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  mae  <- r$mae
  rmse <- r$rmse
  labs <- colnames(mae)
  B    <- nrow(mae)

  keep_mae  <- mae
  keep_rmse <- rmse
  keep_mae [!is.finite(keep_mae)  | keep_mae  > cut_mae ] <- NA_real_
  keep_rmse[!is.finite(keep_rmse) | keep_rmse > cut_rmse] <- NA_real_

  dropped_mae  <- colSums(!is.finite(mae)  | mae  > cut_mae,  na.rm = TRUE)
  dropped_rmse <- colSums(!is.finite(rmse) | rmse > cut_rmse, na.rm = TRUE)

  thr_mae  <- stats::quantile(as.vector(keep_mae),  0.95, na.rm = TRUE)
  thr_rmse <- stats::quantile(as.vector(keep_rmse), 0.95, na.rm = TRUE)

  sum_tbl <- data.frame(
    Model     = labs,
    Mean_MAE  = colMeans(keep_mae,  na.rm = TRUE),
    Mean_RMSE = colMeans(keep_rmse, na.rm = TRUE),
    Max_MAE   = apply(keep_mae,  2, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)),
    Max_RMSE  = apply(keep_rmse, 2, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)),
    Dropped_MAE  = dropped_mae,
    Dropped_RMSE = dropped_rmse,
    DroppedPct_MAE  = dropped_mae  / B,
    DroppedPct_RMSE = dropped_rmse / B,
    Prop_High_MAE  = colMeans(keep_mae  > thr_mae,  na.rm = TRUE),
    Prop_High_RMSE = colMeans(keep_rmse > thr_rmse, na.rm = TRUE)
  )

  out_csv <- file.path(out_dir, paste0("SingleCVB_Summary_", region_label, "_ROBUST.csv"))
  utils::write.csv(sum_tbl, out_csv, row.names = FALSE)

  invisible(list(summary = sum_tbl, mae = keep_mae, rmse = keep_rmse,
                 thr_mae = thr_mae, thr_rmse = thr_rmse))
}

#' Build Table 3.1.2-style per-model summary (MAE/RMSE with bootstrap CIs + success rate)
#'
#' Uses the RDS (source of truth) + robust trimming for reporting. Percentile CIs are
#' computed from the trimmed matrices.
#'
#' @param rds_path path to SingleCVB_Paired_*.rds.
#' @param region_label label used only for output naming.
#' @param out_csv output CSV path.
#' @param cut_mae,cut_rmse robust cutoffs.
#' @param out_dir output directory (optional; used if out_csv is relative).
#' @return data.frame written to out_csv.
make_table_312 <- function(rds_path, region_label, out_csv,
                           cut_mae = 1e3, cut_rmse = 1e3,
                           out_dir = ".") {
  r <- readRDS(rds_path)
  rob <- robust_summarize(r, region_label, cut_mae = cut_mae, cut_rmse = cut_rmse, out_dir = out_dir)
  mae  <- rob$mae
  rmse <- rob$rmse
  labs <- colnames(mae)

  q2.5  <- function(x) stats::quantile(x, 0.025, na.rm = TRUE)
  q97.5 <- function(x) stats::quantile(x, 0.975, na.rm = TRUE)

  ci_tbl <- data.frame(
    Model  = labs,
    MAE_L  = apply(mae,  2, q2.5),
    MAE_U  = apply(mae,  2, q97.5),
    RMSE_L = apply(rmse, 2, q2.5),
    RMSE_U = apply(rmse, 2, q97.5)
  )

  conv <- r$convergence
  conv$SuccessRate <- 100 * conv$SuccessRate_FoldLevel

  out <- merge(rob$summary, ci_tbl, by = "Model", all.x = TRUE)
  out <- merge(out, conv[, c("Model","SuccessfulFits","TotalFitsAttempted","SuccessRate")],
               by = "Model", all.x = TRUE)

  rnd <- function(x, k=2) ifelse(is.finite(x), round(x, k), NA_real_)
  out$`MAE (95% CI)`  <- sprintf("%.2f (%.2f–%.2f)", rnd(out$Mean_MAE),  rnd(out$MAE_L),  rnd(out$MAE_U))
  out$`RMSE (95% CI)` <- sprintf("%.2f (%.2f–%.2f)", rnd(out$Mean_RMSE), rnd(out$RMSE_L), rnd(out$RMSE_U))
  out$`Success rate (%)` <- sprintf("%.2f (%s/%s)", rnd(out$SuccessRate, 2),
                                    format(out$SuccessfulFits, scientific = FALSE),
                                    format(out$TotalFitsAttempted, scientific = FALSE))

  final <- out[, c("Model","MAE (95% CI)","RMSE (95% CI)","Success rate (%)",
                   "Max_MAE","Max_RMSE","DroppedPct_MAE","DroppedPct_RMSE",
                   "Prop_High_MAE","Prop_High_RMSE")]

  out_path <- if (grepl("^(/|[A-Za-z]:)", out_csv)) out_csv else file.path(out_dir, out_csv)
  utils::write.csv(final, out_path, row.names = FALSE)
  invisible(final)
}

#' Build Table 3.1.3-style pairwise comparisons (paired deltas for MAE and RMSE)
#'
#' For each model pair A vs B, computes paired bootstrap deltas on the SAME rows:
#' Δ = A − B; negative Δ means A has lower error. Also returns win% and N_pairs.
#'
#' @param rds_path path to SingleCVB_Paired_*.rds.
#' @param region_label character label (printed in the table).
#' @param out_csv output CSV path.
#' @param cut_mae,cut_rmse robust cutoffs.
#' @param as_percent if TRUE, win fractions are expressed as percent.
#' @param out_dir output directory (optional; used if out_csv is relative).
#' @return data.frame written to out_csv.
make_table_313 <- function(rds_path, region_label, out_csv,
                           cut_mae = 1e3, cut_rmse = 1e3,
                           as_percent = TRUE,
                           out_dir = ".") {
  r <- readRDS(rds_path)
  mae  <- r$mae
  rmse <- r$rmse
  labs <- colnames(mae)

  mae [!is.finite(mae)  | mae  > cut_mae ] <- NA_real_
  rmse[!is.finite(rmse) | rmse > cut_rmse] <- NA_real_

  M <- length(labs)
  out <- list()

  for (i in seq_len(M)) for (j in seq_len(M)) if (i < j) {
    a <- labs[i]; b <- labs[j]

    # paired rows only
    ok_mae  <- is.finite(mae[, a])  & is.finite(mae[, b])
    d_mae   <- mae[ok_mae, a] - mae[ok_mae, b]

    ok_rmse <- is.finite(rmse[, a]) & is.finite(rmse[, b])
    d_rmse  <- rmse[ok_rmse, a] - rmse[ok_rmse, b]

    out[[paste(a, "vs", b)]] <- data.frame(
      Region = region_label,
      Comparison = paste0(sub("^Model\\s*", "M", a), " vs ", sub("^Model\\s*", "M", b)),

      Delta_MAE = if (length(d_mae)) mean(d_mae) else NA_real_,
      CI_L_MAE  = if (length(d_mae)) as.numeric(stats::quantile(d_mae, 0.025)) else NA_real_,
      CI_U_MAE  = if (length(d_mae)) as.numeric(stats::quantile(d_mae, 0.975)) else NA_real_,
      Win_MAE   = if (length(d_mae)) mean(d_mae < 0) else NA_real_,
      N_pairs_MAE = length(d_mae),

      Delta_RMSE = if (length(d_rmse)) mean(d_rmse) else NA_real_,
      CI_L_RMSE  = if (length(d_rmse)) as.numeric(stats::quantile(d_rmse, 0.025)) else NA_real_,
      CI_U_RMSE  = if (length(d_rmse)) as.numeric(stats::quantile(d_rmse, 0.975)) else NA_real_,
      Win_RMSE   = if (length(d_rmse)) mean(d_rmse < 0) else NA_real_,
      N_pairs_RMSE = length(d_rmse),

      stringsAsFactors = FALSE
    )
  }

  tbl <- do.call(rbind, out)

  if (as_percent) {
    tbl$Win_MAE  <- 100 * tbl$Win_MAE
    tbl$Win_RMSE <- 100 * tbl$Win_RMSE
  }

  out_path <- if (grepl("^(/|[A-Za-z]:)", out_csv)) out_csv else file.path(out_dir, out_csv)
  utils::write.csv(tbl, out_path, row.names = FALSE)
  invisible(tbl)
}
