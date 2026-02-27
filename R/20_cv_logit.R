# R/20_cv_logit.R
# Logistic (binomial) nested bootstrap CV with inner AICc selection
# Public template: generic + data-agnostic

#' Safe Hosmer–Lemeshow test with adaptive bin count
#' Returns a list: list(g = <bins used or NA>, hl = <hoslem.test object or NULL>)
.safe_hl <- function(actual, prob, g_max = 10, g_min = 3) {
  uprob <- length(unique(prob))
  if (uprob <= g_min) return(list(g = NA_integer_, hl = NULL))

  g_cap_prob   <- max(g_min, min(g_max, uprob - 1))
  g_cap_sample <- max(g_min, min(g_max, floor(length(actual) / 5)))  # ~5 obs/bin
  g_start <- min(g_cap_prob, g_cap_sample)

  for (g in seq(g_start, g_min, by = -1)) {
    out <- suppressWarnings(
      try(ResourceSelection::hoslem.test(actual, prob, g = g), silent = TRUE)
    )
    if (!inherits(out, "try-error") && !is.null(out$p.value) && is.finite(out$p.value)) {
      return(list(g = g, hl = out))
    }
  }
  list(g = NA_integer_, hl = NULL)
}

#' Simple mode (most frequent value), ignoring NAs
Mode <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_character_)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#' Nested bootstrap CV for logistic regression with inner AICc selection
#'
#' Outer loop: bootstrap resampling
#' Inner loop: leave-one-out within bootstrap sample; select best model by AICc on training fold
#'
#' @param data data.frame containing outcome + predictors
#' @param model_formulas named list of formulas (e.g., list("M1"=y~x1, ...))
#' @param region_label character tag used only for tracking (not required)
#' @param outcome_var name of binary outcome column (0/1 numeric or factor with 2 levels)
#' @param B number of bootstrap iterations
#' @param seed RNG seed
#' @param eps probability clamp to avoid log(0)
#' @return list with metrics, selection frequencies, convergence matrix, and per-bootstrap OOF predictions
nested_bootstrap_CV_logit <- function(
    data,
    model_formulas,
    region_label = "logit_template",
    outcome_var,
    B = 500,
    seed = 123,
    eps = 1e-8
) {
  # packages needed inside function (fail fast with clear messages)
  if (!requireNamespace("MuMIn", quietly = TRUE)) stop("Package 'MuMIn' is required.")
  if (!requireNamespace("pROC", quietly = TRUE)) stop("Package 'pROC' is required.")
  if (!requireNamespace("ResourceSelection", quietly = TRUE)) stop("Package 'ResourceSelection' is required.")

  set.seed(seed)
  n  <- nrow(data)
  mK <- length(model_formulas)
  if (mK < 1) stop("model_formulas must contain at least 1 formula.")

  # ensure named models
  if (is.null(names(model_formulas)) || any(names(model_formulas) == "")) {
    names(model_formulas) <- paste0("Model ", seq_along(model_formulas))
  }

  # storage
  boot_AUC     <- rep(NA_real_, B)
  boot_Brier   <- rep(NA_real_, B)
  boot_LogLoss <- rep(NA_real_, B)
  boot_CalInt  <- rep(NA_real_, B)
  boot_CalSlope<- rep(NA_real_, B)
  boot_HL_stat <- rep(NA_real_, B)
  boot_HL_p    <- rep(NA_real_, B)
  boot_HLg     <- rep(NA_integer_, B)
  boot_SelectedMode <- rep(NA_character_, B)

  conv_mat <- matrix(NA_integer_, nrow = B, ncol = mK,
                     dimnames = list(paste0("B", 1:B), names(model_formulas)))
  all_boot <- vector("list", B)

  # helper: coerce outcome to 0/1 numeric safely
  .coerce_y01 <- function(y) {
    if (is.factor(y)) {
      if (nlevels(y) != 2) stop("Outcome factor must have exactly 2 levels.")
      y <- as.integer(y) - 1L
    } else if (is.logical(y)) {
      y <- as.integer(y)
    } else {
      y <- as.numeric(y)
    }
    if (!all(y %in% c(0, 1, NA))) stop("Outcome must be coded as 0/1 (or a 2-level factor).")
    y
  }

  for (b in 1:B) {
    idx_boot <- sample.int(n, replace = TRUE)
    boot_df  <- data[idx_boot, , drop = FALSE]

    oof_prob   <- rep(NA_real_, n)
    oof_actual <- .coerce_y01(boot_df[[outcome_var]])
    chosen_by_i<- rep(NA_character_, n)

    # LOO within the bootstrap sample
    for (i in 1:n) {
      train_df <- boot_df[-i, , drop = FALSE]
      test_df  <- boot_df[i,  , drop = FALSE]

      # inner selection by AICc on the training fold
      aicc_vec <- rep(Inf, mK)

      for (k in seq_len(mK)) {
        fml <- model_formulas[[k]]
        fit <- try(stats::glm(fml, data = train_df, family = stats::binomial()), silent = TRUE)
        ok  <- !inherits(fit, "try-error") &&
          isTRUE(fit$converged) &&
          all(is.finite(stats::coef(fit)))

        conv_mat[b, k] <- if (ok) 1L else 0L
        if (ok) aicc_vec[k] <- MuMIn::AICc(fit)
      }

      best_k <- which.min(aicc_vec)
      best_f <- model_formulas[[best_k]]
      best_m <- try(stats::glm(best_f, data = train_df, family = stats::binomial()), silent = TRUE)

      if (!inherits(best_m, "try-error") && isTRUE(best_m$converged)) {
        pr <- try(stats::predict(best_m, newdata = test_df, type = "response"), silent = TRUE)
        if (!inherits(pr, "try-error") && length(pr) == 1L && is.finite(pr)) {
          oof_prob[i]    <- as.numeric(pr)
          chosen_by_i[i] <- names(model_formulas)[best_k]
        }
      }
    }

    # metrics on out-of-fold predictions
    keep <- is.finite(oof_prob) & is.finite(oof_actual)

    # AUC only if both classes present
    if (sum(keep) >= 5 && length(unique(oof_actual[keep])) == 2) {
      roc_obj <- pROC::roc(oof_actual[keep], oof_prob[keep], quiet = TRUE)
      boot_AUC[b] <- as.numeric(pROC::auc(roc_obj))
    }

    p <- pmin(pmax(oof_prob, eps), 1 - eps)
    y <- oof_actual
    boot_Brier[b]   <- mean((p - y)^2, na.rm = TRUE)
    boot_LogLoss[b] <- -mean(y * log(p) + (1 - y) * log(1 - p), na.rm = TRUE)

    # calibration intercept/slope (on kept pairs)
    if (sum(keep) >= 5) {
      lp <- stats::qlogis(p[keep])
      cint <- try(stats::coef(stats::glm(y[keep] ~ 1, family = stats::binomial(), offset = lp))[1],
                  silent = TRUE)
      boot_CalInt[b] <- if (!inherits(cint, "try-error")) as.numeric(cint) else NA_real_

      sl <- try(stats::coef(stats::glm(y[keep] ~ lp, family = stats::binomial()))[2],
                silent = TRUE)
      boot_CalSlope[b] <- if (!inherits(sl, "try-error")) as.numeric(sl) else NA_real_
    }

    # Hosmer–Lemeshow (safe)
    hlr <- .safe_hl(y[keep], p[keep], g_max = 10, g_min = 3)
    if (!is.null(hlr$hl)) {
      boot_HLg[b]     <- hlr$g
      boot_HL_stat[b] <- as.numeric(hlr$hl$statistic)
      boot_HL_p[b]    <- as.numeric(hlr$hl$p.value)
    }

    boot_SelectedMode[b] <- Mode(chosen_by_i)
    all_boot[[b]] <- list(oof_prob = oof_prob, oof_actual = oof_actual, chosen_by_i = chosen_by_i)
  }

  # summaries
  ci_tab <- data.frame(
    metric  = c("AUC", "Brier", "LogLoss", "CalIntercept", "CalSlope", "HL_p"),
    mean    = c(mean(boot_AUC, na.rm = TRUE),
                mean(boot_Brier, na.rm = TRUE),
                mean(boot_LogLoss, na.rm = TRUE),
                mean(boot_CalInt, na.rm = TRUE),
                mean(boot_CalSlope, na.rm = TRUE),
                mean(boot_HL_p, na.rm = TRUE)),
    sd      = c(stats::sd(boot_AUC, na.rm = TRUE),
                stats::sd(boot_Brier, na.rm = TRUE),
                stats::sd(boot_LogLoss, na.rm = TRUE),
                stats::sd(boot_CalInt, na.rm = TRUE),
                stats::sd(boot_CalSlope, na.rm = TRUE),
                stats::sd(boot_HL_p, na.rm = TRUE)),
    lower95 = c(stats::quantile(boot_AUC, 0.025, na.rm = TRUE),
                stats::quantile(boot_Brier, 0.025, na.rm = TRUE),
                stats::quantile(boot_LogLoss, 0.025, na.rm = TRUE),
                stats::quantile(boot_CalInt, 0.025, na.rm = TRUE),
                stats::quantile(boot_CalSlope, 0.025, na.rm = TRUE),
                stats::quantile(boot_HL_p, 0.025, na.rm = TRUE)),
    upper95 = c(stats::quantile(boot_AUC, 0.975, na.rm = TRUE),
                stats::quantile(boot_Brier, 0.975, na.rm = TRUE),
                stats::quantile(boot_LogLoss, 0.975, na.rm = TRUE),
                stats::quantile(boot_CalInt, 0.975, na.rm = TRUE),
                stats::quantile(boot_CalSlope, 0.975, na.rm = TRUE),
                stats::quantile(boot_HL_p, 0.975, na.rm = TRUE))
  )

  sel_freq <- sort(table(boot_SelectedMode), decreasing = TRUE)

  list(
    region_label = region_label,
    ci = ci_tab,
    AUC = boot_AUC, Brier = boot_Brier, LogLoss = boot_LogLoss,
    CalInt = boot_CalInt, CalSlope = boot_CalSlope,
    HL_stat = boot_HL_stat, HL_p = boot_HL_p, HL_g = boot_HLg,
    SelectedMode = boot_SelectedMode,
    SelectionFreq = sel_freq,
    Convergence = conv_mat,
    PerBootstrap = all_boot
  )
}
