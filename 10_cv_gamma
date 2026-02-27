# R/10_cv_gamma.R
# Nested bootstrap LOOCV for Gamma(log) with inner AICc model selection
# - fold-wise scaling (numeric predictors only) to prevent leakage
# - convergence logging (inner + outer refit)
# - guardrails for invalid predictions
# - returns a results list; saving happens outside this function

suppressPackageStartupMessages({
  library(MuMIn)   # AICc
  library(Metrics) # mae, rmse
})

# --- helpers ---------------------------------------------------------

Mode <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_character_)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

get_numeric_terms <- function(df, formula) {
  trm <- attr(terms(formula), "term.labels")
  if (!length(trm)) return(character(0))
  trm[vapply(trm, function(v) is.numeric(df[[v]]), logical(1))]
}

scale_train <- function(df, vars) {
  if (!length(vars)) return(list(df = df, center = NULL, scale = NULL))
  cen <- vapply(vars, \(v) mean(df[[v]], na.rm = TRUE), numeric(1))
  sc  <- vapply(vars, \(v) sd(df[[v]],   na.rm = TRUE), numeric(1))
  sc[!is.finite(sc) | sc == 0] <- 1
  d2 <- df
  for (v in vars) d2[[v]] <- (df[[v]] - cen[v]) / sc[v]
  list(df = d2, center = cen, scale = sc)
}

scale_apply_row <- function(row, center, scale) {
  if (is.null(center)) return(row)
  for (v in names(center)) row[[v]] <- (row[[v]] - center[v]) / scale[v]
  row
}

# --- main engine ------------------------------------------------------

nested_bootstrap_cv_gamma <- function(
    data,
    model_formulas,
    outcome_var,
    B = 500,
    seed = 123,
    scale_predictors = TRUE,
    guardrail = c("na", "clip"),
    min_y = 1,
    max_y = 900,
    high_error_threshold = 300
) {
  guardrail <- match.arg(guardrail)
  set.seed(seed)

  stopifnot(is.data.frame(data))
  stopifnot(outcome_var %in% names(data))
  stopifnot(length(model_formulas) >= 1)
  if (is.null(names(model_formulas)) || any(names(model_formulas) == "")) {
    stop("model_formulas must be a *named* list (e.g., list('Model 1' = y ~ x1, ...)).")
  }

  n <- nrow(data)
  M <- length(model_formulas)

  bootstrap_mae  <- rep(NA_real_, B)
  bootstrap_rmse <- rep(NA_real_, B)
  bootstrap_selected_model <- rep(NA_character_, B)

  # proportion of inner fits that converged: bootstrap x model
  convergence_matrix <- matrix(NA_real_, nrow = B, ncol = M,
                               dimnames = list(paste0("b", seq_len(B)), names(model_formulas)))

  # detailed inner fold convergence: bootstrap x fold x model (0/1)
  conv_detail <- array(NA_integer_, dim = c(B, n, M),
                       dimnames = list(paste0("b", seq_len(B)),
                                       paste0("i", seq_len(n)),
                                       names(model_formulas)))

  # outer refit convergence: bootstrap x fold (0/1) + chosen model label
  outer_refit_ok   <- matrix(NA_integer_,   nrow = B, ncol = n)
  outer_best_model <- matrix(NA_character_, nrow = B, ncol = n)

  failed_log <- list()
  per_bootstrap_results <- vector("list", B)

  for (b in seq_len(B)) {

    idx_boot <- sample.int(n, replace = TRUE)
    boot_data <- data[idx_boot, , drop = FALSE]

    outer_pred <- rep(NA_real_, n)
    outer_y    <- boot_data[[outcome_var]]
    outer_chosen <- rep(NA_character_, n)

    success_counts <- numeric(M)

    for (i in seq_len(n)) {
      outer_train <- boot_data[-i, , drop = FALSE]
      outer_test  <- boot_data[i,  , drop = FALSE]

      inner_aicc <- rep(Inf, M)

      for (m in seq_along(model_formulas)) {
        fml <- model_formulas[[m]]
        tr_df <- outer_train

        if (scale_predictors) {
          vars <- get_numeric_terms(tr_df, fml)
          s <- scale_train(tr_df, vars)
          tr_df <- s$df
        }

        fit <- try(glm(fml, data = tr_df, family = Gamma(link = "log")), silent = TRUE)
        ok  <- (!inherits(fit, "try-error")) && isTRUE(fit$converged)

        conv_detail[b, i, m] <- as.integer(ok)

        if (ok) {
          inner_aicc[m] <- MuMIn::AICc(fit)
          success_counts[m] <- success_counts[m] + 1
        } else {
          inner_aicc[m] <- Inf
        }
      }

      best_idx <- which.min(inner_aicc)
      best_fml <- model_formulas[[best_idx]]
      best_lbl <- names(model_formulas)[best_idx]

      outer_chosen[i] <- best_lbl
      outer_best_model[b, i] <- best_lbl

      # refit best on outer train (with fold-wise scaling), predict outer test
      tr_df <- outer_train
      te_df <- outer_test

      if (scale_predictors) {
        vars <- get_numeric_terms(tr_df, best_fml)
        s <- scale_train(tr_df, vars)
        tr_df <- s$df
        te_df <- scale_apply_row(te_df, s$center, s$scale)
      }

      best_fit <- try(glm(best_fml, data = tr_df, family = Gamma(link = "log")), silent = TRUE)
      ok_outer <- (!inherits(best_fit, "try-error")) && isTRUE(best_fit$converged)
      outer_refit_ok[b, i] <- as.integer(ok_outer)

      if (ok_outer) {
        pr <- try(predict(best_fit, newdata = te_df, type = "response"), silent = TRUE)
        if (!inherits(pr, "try-error")) {
          p <- as.numeric(pr)
          bad <- (!is.finite(p)) || p <= 0 || p > max_y
          if (bad) {
            if (guardrail == "clip") {
              p <- min(max(p, min_y), max_y)
              outer_pred[i] <- p
            } else {
              outer_pred[i] <- NA_real_
            }
          } else {
            outer_pred[i] <- p
          }
        }
      } else {
        failed_log[[length(failed_log) + 1]] <- list(
          bootstrap = b,
          fold = i,
          model = best_lbl,
          formula = deparse(best_fml)
        )
      }
    } # end outer LOOCV

    keep <- is.finite(outer_pred) & is.finite(outer_y)

    bootstrap_mae[b]  <- if (any(keep)) Metrics::mae (outer_y[keep], outer_pred[keep]) else NA_real_
    bootstrap_rmse[b] <- if (any(keep)) Metrics::rmse(outer_y[keep], outer_pred[keep]) else NA_real_

    bootstrap_selected_model[b] <- Mode(outer_chosen)

    # divide by n outer folds (proportion of successful inner fits)
    convergence_matrix[b, ] <- success_counts / n

    per_bootstrap_results[[b]] <- list(
      mae = bootstrap_mae[b],
      rmse = bootstrap_rmse[b],
      selected = bootstrap_selected_model[b],
      all_models = outer_chosen,
      actual = outer_y,
      predicted = outer_pred
    )
  } # end B

  ci_summary <- data.frame(
    metric  = c("MAE", "RMSE"),
    mean    = c(mean(bootstrap_mae,  na.rm = TRUE), mean(bootstrap_rmse, na.rm = TRUE)),
    sd      = c(sd(bootstrap_mae,    na.rm = TRUE), sd(bootstrap_rmse,   na.rm = TRUE)),
    lower95 = c(quantile(bootstrap_mae,  0.025, na.rm = TRUE),
                quantile(bootstrap_rmse, 0.025, na.rm = TRUE)),
    upper95 = c(quantile(bootstrap_mae,  0.975, na.rm = TRUE),
                quantile(bootstrap_rmse, 0.975, na.rm = TRUE))
  )

  robust_summary <- data.frame(
    metric = c("MAE", "RMSE"),
    median = c(median(bootstrap_mae,  na.rm = TRUE), median(bootstrap_rmse, na.rm = TRUE)),
    IQR_L  = c(quantile(bootstrap_mae,  0.25, na.rm = TRUE),
               quantile(bootstrap_rmse, 0.25, na.rm = TRUE)),
    IQR_H  = c(quantile(bootstrap_mae,  0.75, na.rm = TRUE),
               quantile(bootstrap_rmse, 0.75, na.rm = TRUE)),
    p95    = c(quantile(bootstrap_mae,  0.95, na.rm = TRUE),
               quantile(bootstrap_rmse, 0.95, na.rm = TRUE))
  )

  prop_high_mae <- mean(bootstrap_mae > high_error_threshold, na.rm = TRUE)
  top_model_per_iter <- vapply(per_bootstrap_results, `[[`, character(1), "selected")

  list(
    settings = list(
      B = B, seed = seed, outcome_var = outcome_var,
      scale_predictors = scale_predictors,
      guardrail = guardrail, min_y = min_y, max_y = max_y,
      high_error_threshold = high_error_threshold
    ),
    mae = bootstrap_mae,
    rmse = bootstrap_rmse,
    ci_summary = ci_summary,
    robust_summary = robust_summary,
    prop_high_mae = prop_high_mae,
    model_frequencies = table(bootstrap_selected_model),
    failed_log = failed_log,
    convergence_matrix = convergence_matrix,
    conv_detail = conv_detail,
    outer_refit_ok = outer_refit_ok,
    outer_best_model = outer_best_model,
    per_bootstrap_results = per_bootstrap_results,
    top_model_per_iter = top_model_per_iter
  )
}
