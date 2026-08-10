set.seed(20260309)

suppressPackageStartupMessages({
  library(survey)
  library(ggplot2)
})

weighted_cov_matrix <- function(mat, w) {
  w <- w / sum(w)
  mu <- colSums(mat * w)
  centered <- sweep(mat, 2, mu, "-")
  zw <- centered * sqrt(w)
  crossprod(zw)
}

simulate_population <- function(N_pop = 100000, p_vars = 25, target_rr = 0.30) {
  z_rr <- rnorm(N_pop)
  z_aux <- rnorm(N_pop)

  X <- matrix(0L, nrow = N_pop, ncol = p_vars)
  for (j in seq_len(p_vars)) {
    if (j <= 8) {
      eta <- runif(1, 0.9, 1.4) * z_rr + 0.2 * z_aux + rnorm(N_pop, sd = 0.70)
    } else if (j <= 16) {
      eta <- runif(1, 0.4, 0.8) * z_rr + 0.4 * z_aux + rnorm(N_pop, sd = 1.00)
    } else {
      eta <- runif(1, 0.1, 0.4) * z_rr + 0.6 * z_aux + rnorm(N_pop, sd = 1.10)
    }
    shift <- runif(1, -0.5, 0.5)
    X[, j] <- rbinom(N_pop, 1, plogis(shift + eta))
  }
  colnames(X) <- paste0("X", seq_len(p_vars))

  rr_coefs <- c(seq(0.9, 0.35, length.out = 8),
                0.30, 0.25, -0.20, 0.18,
                rep(0.08, p_vars - 12))
  lin_rr <- as.vector(X %*% rr_coefs + 0.30 * z_rr)
  rr_intercept <- uniroot(function(a) mean(plogis(a + lin_rr)) - target_rr,
                          interval = c(-10, 4))$root
  rho <- plogis(rr_intercept + lin_rr)

  y_coefs <- seq(0.9, 0.08, length.out = p_vars)
  interaction_pairs <- matrix(
    c(
      1, 2,
      1, 3,
      2, 4,
      3, 5,
      4, 6,
      5, 7,
      6, 8,
      1, 8,
      2, 7,
      3, 6,
      4, 5,
      9, 10,
      11, 12,
      13, 14
    ),
    byrow = TRUE, ncol = 2
  )
  interaction_coefs <- c(1.2, -1.0, 0.9, -0.8, 0.7, -0.6, 0.55,
                         -0.45, 0.4, -0.35, 0.3, 0.25, -0.2, 0.18)
  interaction_signal <- rowSums(
    vapply(seq_len(nrow(interaction_pairs)), function(ii) {
      p1 <- interaction_pairs[ii, 1]
      p2 <- interaction_pairs[ii, 2]
      interaction_coefs[ii] * (X[, p1] * X[, p2])
    }, numeric(N_pop))
  )

  y <- 0.9 * z_rr + 0.35 * z_aux + as.vector(X %*% y_coefs) +
    interaction_signal + rnorm(N_pop, sd = 1.1)

  list(
    X = X,
    rho = rho,
    y = y,
    z_rr = z_rr,
    interaction_pairs = interaction_pairs,
    interaction_signal = interaction_signal
  )
}

safe_calibrate <- function(...) {
  failed <- FALSE
  fit <- tryCatch(
    withCallingHandlers(
      calibrate(...),
      warning = function(w) {
        if (grepl("Failed to converge", conditionMessage(w), fixed = TRUE)) {
          failed <<- TRUE
        }
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) NULL
  )
  if (is.null(fit) || failed) {
    return(NULL)
  }
  fit
}

run_simulation <- function(
  reps = 100,
  sample_n = 7000,
  K_values = 2:25,
  N_pop = 100000,
  p_vars = 25,
  target_rr = 0.30,
  eval_n = 12000
) {
  pop <- simulate_population(N_pop = N_pop, p_vars = p_vars, target_rr = target_rr)
  X <- pop$X
  rho <- pop$rho
  y <- pop$y
  interaction_pairs <- pop$interaction_pairs
  interaction_signal <- pop$interaction_signal

  true_means <- colMeans(X)
  true_y_mean <- mean(y)
  true_cov_X <- lapply(K_values, function(k) cov(X[, seq_len(k), drop = FALSE]))
  names(true_cov_X) <- as.character(K_values)

  eval_idx <- sample.int(N_pop, eval_n, replace = FALSE)
  eval_df <- as.data.frame(X[eval_idx, , drop = FALSE])
  eval_y <- y[eval_idx]

  corr_x_rho <- vapply(seq_len(p_vars), function(j) cor(X[, j], rho), numeric(1))
  cat("Population diagnostics:\n")
  cat(sprintf("  Mean response rate: %.3f\n", mean(rho)))
  cat(sprintf("  Mean |cor(X_j, rho)| for j=1..8: %.3f\n",
              mean(abs(corr_x_rho[1:8]))))
  cat(sprintf("  Mean |cor(X_j, rho)| for j=9..25: %.3f\n",
              mean(abs(corr_x_rho[9:p_vars]))))
  cat(sprintf("  Cor(y, interaction_signal): %.3f\n",
              cor(y, interaction_signal)))

  out <- vector("list", length = reps * length(K_values) * 2)
  row_id <- 1L

  for (r in seq_len(reps)) {
    sample_idx <- sample.int(N_pop, sample_n, replace = FALSE)
    responded <- rbinom(sample_n, 1, rho[sample_idx]) == 1
    rr_obs <- mean(responded)
    resp_idx <- sample_idx[responded]
    n_resp <- length(resp_idx)

    if (n_resp < 400) {
      next
    }

    df <- as.data.frame(X[resp_idx, , drop = FALSE])
    df$y <- y[resp_idx]
    base_design <- svydesign(ids = ~1, data = df, weights = ~1)

    for (k in K_values) {
      vars <- paste0("X", seq_len(k))
      form <- as.formula(paste("~", paste(vars, collapse = " + ")))
      pop_totals <- c("(Intercept)" = n_resp, setNames(n_resp * true_means[seq_len(k)], vars))

      std <- safe_calibrate(
        base_design,
        formula = form,
        population = pop_totals,
        calfun = "raking",
        maxit = 350,
        epsilon = 1e-8
      )

      con <- safe_calibrate(
        base_design,
        formula = form,
        population = pop_totals,
        calfun = "raking",
        bounds = c(rr_obs, Inf),
        maxit = 350,
        epsilon = 1e-8
      )

      if (is.null(std) || is.null(con)) {
        next
      }

      for (method in c("standard", "constrained")) {
        design_obj <- if (method == "standard") std else con
        w <- as.numeric(weights(design_obj))
        mat <- as.matrix(df[, vars, drop = FALSE])

        cov_hat <- weighted_cov_matrix(mat, w)
        truth <- true_cov_X[[as.character(k)]]
        off_diag <- upper.tri(cov_hat)
        cov_mae_x <- mean(abs(cov_hat[off_diag] - truth[off_diag]))

        w_norm <- w / sum(w)
        margin_hat <- colSums(mat * w_norm)
        margin_max_abs_err <- max(abs(margin_hat - true_means[seq_len(k)]))

        y_mean_hat <- sum(w_norm * df$y)
        y_mean_abs_err <- abs(y_mean_hat - true_y_mean)

        fit_df <- df[, c("y", vars), drop = FALSE]
        active_interactions <- interaction_pairs[
          interaction_pairs[, 1] <= k & interaction_pairs[, 2] <= k,
          ,
          drop = FALSE
        ]
        interaction_terms <- character(0)
        if (nrow(active_interactions) > 0) {
          interaction_terms <- apply(active_interactions, 1, function(pair) {
            paste0("X", pair[1], ":X", pair[2])
          })
        }
        rhs_terms <- c(vars, interaction_terms)
        model_formula <- as.formula(paste("y ~", paste(rhs_terms, collapse = " + ")))

        fit <- tryCatch(
          lm(model_formula,
             data = fit_df, weights = w),
          error = function(e) NULL
        )

        y_pred_mae <- NA_real_
        if (!is.null(fit)) {
          pred <- tryCatch(
            predict(fit, newdata = eval_df[, vars, drop = FALSE]),
            error = function(e) rep(NA_real_, length(eval_y))
          )
          if (!all(is.na(pred))) {
            y_pred_mae <- mean(abs(pred - eval_y))
          }
        }

        out[[row_id]] <- data.frame(
          rep = r,
          K = k,
          method = method,
          rr_obs = rr_obs,
          weight_var = var(w),
          min_weight = min(w),
          max_weight = max(w),
          margin_max_abs_err = margin_max_abs_err,
          cov_mae_x = cov_mae_x,
          y_mean_abs_err = y_mean_abs_err,
          y_pred_mae = y_pred_mae
        )
        row_id <- row_id + 1L
      }
    }
  }

  res <- do.call(rbind, out)
  res <- res[complete.cases(res), ]
  res
}

summarize_results <- function(res) {
  metrics <- c(
    "weight_var", "cov_mae_x", "y_mean_abs_err", "y_pred_mae",
    "margin_max_abs_err", "min_weight", "max_weight"
  )

  all_summaries <- lapply(metrics, function(metric) {
    s <- aggregate(res[[metric]], list(method = res$method, K = res$K), function(x) {
      c(mean = mean(x), median = median(x), sd = sd(x))
    })
    s <- do.call(data.frame, s)
    names(s)[3:5] <- paste0(metric, c("_mean", "_median", "_sd"))
    s
  })

  merged <- Reduce(function(a, b) merge(a, b, by = c("method", "K")), all_summaries)
  merged[order(merged$K, merged$method), ]
}

plot_metric <- function(res, yvar, ylab, file_name) {
  p <- ggplot(res, aes(x = K, y = .data[[yvar]], color = method)) +
    stat_summary(fun = mean, geom = "line", linewidth = 0.9) +
    stat_summary(fun = mean, geom = "point", size = 1.7) +
    labs(x = "Number of weighting variables (K)", y = ylab, color = "") +
    theme_minimal(base_size = 12)
  ggsave(file_name, p, width = 8.2, height = 5.0, dpi = 160)
}

results <- run_simulation()
summary_table <- summarize_results(results)

write.csv(results, "constraint_simulation_long.csv", row.names = FALSE)
write.csv(summary_table, "constraint_simulation_summary.csv", row.names = FALSE)

plot_metric(results, "weight_var", "Variance of calibrated weights", "constraint_sim_weight_var.png")
plot_metric(results, "cov_mae_x", "MAE of pairwise covariances among weighting vars", "constraint_sim_cov_mae_x.png")
plot_metric(results, "y_mean_abs_err", "Absolute error of population mean(y)", "constraint_sim_y_mean_abs_err.png")
plot_metric(results, "y_pred_mae", "Prediction MAE for y (held-out population)", "constraint_sim_y_pred_mae.png")

cat("Completed simulation.\n")
cat(sprintf("Rows in long results: %s\n", nrow(results)))
cat("Summary (means) for key metrics:\n")
print(
  summary_table[, c("method", "K", "weight_var_mean", "cov_mae_x_mean", "y_mean_abs_err_mean", "y_pred_mae_mean")],
  row.names = FALSE
)
