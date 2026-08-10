set.seed(20260310)

suppressPackageStartupMessages({
  library(survey)
  library(ggplot2)
})

weighted_mean <- function(x, w) {
  wn <- w / sum(w)
  sum(wn * x)
}

weighted_cov_binary <- function(x1, x2, w) {
  wn <- w / sum(w)
  m1 <- sum(wn * x1)
  m2 <- sum(wn * x2)
  sum(wn * (x1 - m1) * (x2 - m2))
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

simulate_population <- function(N_pop = 250000) {
  # Joint distribution for (x1, x2) with substantial mass in the interaction cell.
  levels_cell <- c("00", "10", "01", "11")
  probs_cell <- c(0.30, 0.22, 0.22, 0.26)
  cell <- sample(levels_cell, size = N_pop, replace = TRUE, prob = probs_cell)

  x1 <- as.integer(substr(cell, 1, 1))
  x2 <- as.integer(substr(cell, 2, 2))

  # Outcome: positive main effects, strong negative interaction.
  y <- 0.6 +
    1.35 * x1 +
    1.35 * x2 -
    3.80 * (x1 * x2) +
    rnorm(N_pop, sd = 1.15)

  # Response propensities: interaction cell (11) is heavily undersampled.
  rr_map <- c("00" = 0.46, "10" = 0.48, "01" = 0.54, "11" = 0.01)
  rho <- unname(rr_map[cell])

  list(
    x1 = x1,
    x2 = x2,
    y = y,
    cell = factor(cell, levels = levels_cell),
    rho = rho,
    rr_map = rr_map
  )
}

run_simulation <- function(
  reps = 400,
  sample_n = 5000,
  N_pop = 250000,
  eval_n = 20000
) {
  pop <- simulate_population(N_pop = N_pop)

  true_y_mean <- mean(pop$y)
  true_x1_mean <- mean(pop$x1)
  true_x2_mean <- mean(pop$x2)
  true_cov <- cov(pop$x1, pop$x2)
  true_cell <- prop.table(table(pop$cell))
  levels_cell <- names(true_cell)

  eval_idx <- sample.int(N_pop, eval_n, replace = FALSE)
  eval_df <- data.frame(
    x1 = pop$x1[eval_idx],
    x2 = pop$x2[eval_idx],
    y = pop$y[eval_idx]
  )

  cat("Population diagnostics:\n")
  print(round(true_cell, 4))
  cat(sprintf("Target E[y]: %.4f\n", true_y_mean))
  cat(sprintf("Target response rate (mean rho): %.4f\n", mean(pop$rho)))
  cat("Response rates by cell:\n")
  print(pop$rr_map)

  out <- vector("list", length = reps * 5)
  row_id <- 1L

  counts <- list(
    attempted = 0L,
    kept = 0L,
    std_fail = 0L,
    con_fail = 0L,
    oracle_joint_fail = 0L
  )

  form <- ~x1 + x2

  for (r in seq_len(reps)) {
    counts$attempted <- counts$attempted + 1L

    s_idx <- sample.int(N_pop, sample_n, replace = FALSE)
    responded <- rbinom(sample_n, 1, pop$rho[s_idx]) == 1
    rr_obs <- mean(responded)
    ridx <- s_idx[responded]

    if (length(ridx) < 250) {
      next
    }

    df <- data.frame(
      x1 = pop$x1[ridx],
      x2 = pop$x2[ridx],
      y = pop$y[ridx],
      cell = pop$cell[ridx]
    )

    obs_cell <- prop.table(table(factor(df$cell, levels = levels_cell)))
    undersample_ratio_11 <- as.numeric(obs_cell["11"] / true_cell["11"])

    des <- svydesign(ids = ~1, data = df, weights = ~1)
    pop_totals <- c(
      "(Intercept)" = nrow(df),
      x1 = nrow(df) * true_x1_mean,
      x2 = nrow(df) * true_x2_mean
    )

    std <- safe_calibrate(
      des,
      formula = form,
      population = pop_totals,
      calfun = "raking",
      maxit = 500,
      epsilon = 1e-9
    )
    con <- safe_calibrate(
      des,
      formula = form,
      population = pop_totals,
      calfun = "raking",
      bounds = c(rr_obs, Inf),
      maxit = 500,
      epsilon = 1e-9
    )

    if (is.null(std)) {
      counts$std_fail <- counts$std_fail + 1L
    }
    if (is.null(con)) {
      counts$con_fail <- counts$con_fail + 1L
    }

    method_weights <- list(
      unweighted = rep(1, nrow(df)),
      oracle_ipw = {
        w <- 1 / pop$rho[ridx]
        w / mean(w)
      }
    )

    if (!any(obs_cell == 0)) {
      w_joint <- as.numeric(true_cell[as.character(df$cell)] / obs_cell[as.character(df$cell)])
      w_joint <- w_joint / mean(w_joint)
      method_weights$oracle_joint <- w_joint
    } else {
      counts$oracle_joint_fail <- counts$oracle_joint_fail + 1L
    }

    if (!is.null(std)) {
      method_weights$standard <- as.numeric(weights(std))
    }
    if (!is.null(con)) {
      method_weights$constrained <- as.numeric(weights(con))
    }

    for (method in names(method_weights)) {
      w <- method_weights[[method]]
      wn <- w / sum(w)

      y_hat <- sum(wn * df$y)
      cov_hat <- weighted_cov_binary(df$x1, df$x2, w)
      p11_hat <- sum(wn * as.numeric(df$cell == "11"))

      fit <- tryCatch(
        lm(y ~ x1 + x2 + x1:x2, data = df, weights = w),
        error = function(e) NULL
      )
      y_pred_mae <- if (is.null(fit)) {
        NA_real_
      } else {
        pred <- predict(fit, newdata = eval_df)
        mean(abs(pred - eval_df$y))
      }

      out[[row_id]] <- data.frame(
        rep = r,
        method = method,
        rr_obs = rr_obs,
        n_resp = nrow(df),
        undersample_ratio_11 = undersample_ratio_11,
        y_mean_bias = y_hat - true_y_mean,
        y_mean_abs_err = abs(y_hat - true_y_mean),
        y_pred_mae = y_pred_mae,
        cov_bias = cov_hat - true_cov,
        cov_abs_err = abs(cov_hat - true_cov),
        p11_bias = p11_hat - true_cell["11"],
        p11_abs_err = abs(p11_hat - true_cell["11"]),
        weight_var = var(w),
        min_weight = min(w),
        max_weight = max(w)
      )
      row_id <- row_id + 1L
    }

    counts$kept <- counts$kept + 1L
  }

  res <- do.call(rbind, out)
  res <- res[complete.cases(res), ]

  list(results = res, counts = counts)
}

summarize_results <- function(res) {
  aggregate(
    cbind(
      rr_obs, n_resp, undersample_ratio_11, y_mean_bias, y_mean_abs_err, y_pred_mae,
      cov_bias, cov_abs_err, p11_bias, p11_abs_err,
      weight_var, min_weight, max_weight
    ) ~ method,
    data = res,
    FUN = mean
  )
}

compare_methods <- function(res, method_a, method_b) {
  a <- subset(res, method == method_a)
  b <- subset(res, method == method_b)
  paired <- merge(a, b, by = "rep", suffixes = c("_a", "_b"))
  data.frame(
    method_a = method_a,
    method_b = method_b,
    pairs = nrow(paired),
    y_bias_diff_b_minus_a = mean(paired$y_mean_bias_b - paired$y_mean_bias_a, na.rm = TRUE),
    y_abs_err_diff_b_minus_a = mean(paired$y_mean_abs_err_b - paired$y_mean_abs_err_a, na.rm = TRUE),
    y_pred_mae_diff_b_minus_a = mean(paired$y_pred_mae_b - paired$y_pred_mae_a, na.rm = TRUE),
    cov_abs_err_diff_b_minus_a = mean(paired$cov_abs_err_b - paired$cov_abs_err_a, na.rm = TRUE),
    p11_abs_err_diff_b_minus_a = mean(paired$p11_abs_err_b - paired$p11_abs_err_a, na.rm = TRUE),
    weight_var_pct_diff_b_minus_a = mean(
      100 * (paired$weight_var_b - paired$weight_var_a) / paired$weight_var_a,
      na.rm = TRUE
    ),
    share_y_abs_err_improved_b_vs_a = mean(
      paired$y_mean_abs_err_b < paired$y_mean_abs_err_a,
      na.rm = TRUE
    )
  )
}

plot_metric <- function(res, metric, ylab, path) {
  p <- ggplot(res, aes(x = method, y = .data[[metric]], fill = method)) +
    geom_boxplot(width = 0.62, outlier.alpha = 0.25) +
    labs(x = "", y = ylab) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
  ggsave(path, p, width = 7.4, height = 4.8, dpi = 160)
}

sim <- run_simulation()
results <- sim$results
summary_table <- summarize_results(results)
delta <- rbind(
  compare_methods(results, "unweighted", "standard"),
  compare_methods(results, "standard", "constrained"),
  compare_methods(results, "standard", "oracle_joint"),
  compare_methods(results, "constrained", "oracle_joint"),
  compare_methods(results, "standard", "oracle_ipw"),
  compare_methods(results, "constrained", "oracle_ipw")
)

write.csv(results, "negative_interaction_undersampled_long.csv", row.names = FALSE)
write.csv(summary_table, "negative_interaction_undersampled_summary.csv", row.names = FALSE)
write.csv(delta, "negative_interaction_undersampled_delta.csv", row.names = FALSE)

plot_metric(
  results,
  "y_mean_abs_err",
  "Absolute error of mean(y)",
  "negative_interaction_undersampled_abs_err_box.png"
)
plot_metric(
  results,
  "y_mean_bias",
  "Bias of mean(y)",
  "negative_interaction_undersampled_bias_box.png"
)
plot_metric(
  results,
  "p11_abs_err",
  "Absolute error of P(X1=1, X2=1)",
  "negative_interaction_undersampled_p11_err_box.png"
)
plot_metric(
  results,
  "weight_var",
  "Weight variance",
  "negative_interaction_undersampled_weight_var_box.png"
)

cat("Completed negative-interaction undersampled-cell simulation.\n")
cat("Replication counts:\n")
print(sim$counts)
cat("\nMethod means:\n")
print(summary_table, row.names = FALSE)
cat("\nPairwise deltas (method_b - method_a):\n")
print(delta, row.names = FALSE)
