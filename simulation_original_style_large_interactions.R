set.seed(20260310)

suppressPackageStartupMessages({
  library(survey)
  library(ggplot2)
})

weighted_mean <- function(y, w) {
  w <- w / sum(w)
  sum(w * y)
}

simulate_population <- function(N_pop = 200000, K = 8, target_rr = 0.28) {
  z <- rnorm(N_pop)
  X <- matrix(0L, nrow = N_pop, ncol = K)
  for (j in seq_len(K)) {
    # Shared latent z induces substantial covariance among predictors.
    X[, j] <- rbinom(
      N_pop,
      1,
      plogis(runif(1, 0.6, 1.1) * z + rnorm(N_pop, sd = 0.9) + runif(1, -0.4, 0.4))
    )
  }
  colnames(X) <- paste0("x", seq_len(K))

  lin_rr <- -2 +
    0.4 * X[, 1] + 0.4 * X[, 2] + 0.35 * X[, 3] + 0.3 * X[, 4] +
    1.7 * (X[, 1] * X[, 2]) - 1.5 * (X[, 3] * X[, 4]) +
    1.3 * (X[, 5] * X[, 6]) - 1.1 * (X[, 7] * X[, 8]) +
    1.1 * (X[, 1] * X[, 3] * X[, 5]) -
    0.9 * (X[, 2] * X[, 4] * X[, 6])

  intercept <- uniroot(function(a) mean(plogis(a + lin_rr)) - target_rr, c(-8, 4))$root
  rho <- plogis(intercept + lin_rr)

  # Large interaction effects in y (closer to the original examples, but stronger).
  y <- 0.3 +
    0.2 * rowSums(X[, 1:4, drop = FALSE]) +
    2.2 * (X[, 1] * X[, 2]) - 1.8 * (X[, 3] * X[, 4]) +
    2.0 * (X[, 5] * X[, 6]) - 1.7 * (X[, 7] * X[, 8]) +
    1.5 * (X[, 1] * X[, 3] * X[, 5]) -
    1.3 * (X[, 2] * X[, 4] * X[, 6]) +
    rnorm(N_pop, sd = 1.5)

  list(X = X, y = y, rho = rho)
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
  reps = 300,
  n_selected = 3000,
  K = 8,
  N_pop = 200000,
  target_rr = 0.28
) {
  pop <- simulate_population(N_pop = N_pop, K = K, target_rr = target_rr)
  X <- pop$X
  y <- pop$y
  rho <- pop$rho

  true_y_mean <- mean(y)
  true_margins <- colMeans(X)

  out <- vector("list", length = reps * 4)
  row_id <- 1L

  vars <- paste0("x", seq_len(K))
  form <- as.formula(paste("~", paste(vars, collapse = " + ")))

  for (r in seq_len(reps)) {
    s_idx <- sample.int(nrow(X), n_selected, replace = FALSE)
    responded <- rbinom(n_selected, 1, rho[s_idx]) == 1
    rr_obs <- mean(responded)
    ridx <- s_idx[responded]

    if (length(ridx) < 450) {
      next
    }

    df <- as.data.frame(X[ridx, , drop = FALSE])
    names(df) <- vars
    df$y <- y[ridx]

    des <- svydesign(ids = ~1, data = df, weights = ~1)
    pop_totals <- c("(Intercept)" = nrow(df), setNames(nrow(df) * true_margins, vars))

    std <- safe_calibrate(
      des,
      formula = form,
      population = pop_totals,
      calfun = "raking",
      maxit = 300,
      epsilon = 1e-8
    )
    con <- safe_calibrate(
      des,
      formula = form,
      population = pop_totals,
      calfun = "raking",
      bounds = c(rr_obs, Inf),
      maxit = 300,
      epsilon = 1e-8
    )
    if (is.null(std) || is.null(con)) {
      next
    }

    weights_list <- list(
      unweighted = rep(1, nrow(df)),
      standard = as.numeric(weights(std)),
      constrained = as.numeric(weights(con)),
      oracle_ipw = {
        w <- (1 / rho[ridx])
        w / mean(w)
      }
    )

    for (method in names(weights_list)) {
      w <- weights_list[[method]]
      y_hat <- weighted_mean(df$y, w)

      out[[row_id]] <- data.frame(
        rep = r,
        method = method,
        rr_obs = rr_obs,
        y_mean_bias = y_hat - true_y_mean,
        y_mean_abs_err = abs(y_hat - true_y_mean),
        weight_var = var(w),
        min_weight = min(w),
        max_weight = max(w)
      )
      row_id <- row_id + 1L
    }
  }

  res <- do.call(rbind, out)
  res <- res[complete.cases(res), ]
  res
}

summarize_results <- function(res) {
  aggregate(
    cbind(rr_obs, y_mean_bias, y_mean_abs_err, weight_var, min_weight, max_weight) ~ method,
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
    mean_rr = mean(paired$rr_obs_a),
    bias_diff_b_minus_a = mean(paired$y_mean_bias_b - paired$y_mean_bias_a),
    abs_err_diff_b_minus_a = mean(paired$y_mean_abs_err_b - paired$y_mean_abs_err_a),
    weight_var_pct_diff_b_minus_a = mean(100 * (paired$weight_var_b - paired$weight_var_a) / paired$weight_var_a),
    share_abs_err_improved_b_vs_a = mean(paired$y_mean_abs_err_b < paired$y_mean_abs_err_a)
  )
}

plot_metric <- function(res, metric, ylab, path) {
  p <- ggplot(res, aes(x = method, y = .data[[metric]], fill = method)) +
    geom_boxplot(width = 0.6, outlier.alpha = 0.25) +
    labs(x = "", y = ylab) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
  ggsave(path, p, width = 7.4, height = 4.8, dpi = 160)
}

results <- run_simulation()
summary_table <- summarize_results(results)
delta <- rbind(
  compare_methods(results, "unweighted", "standard"),
  compare_methods(results, "standard", "constrained"),
  compare_methods(results, "standard", "oracle_ipw"),
  compare_methods(results, "constrained", "oracle_ipw")
)

write.csv(results, "original_style_interactions_long.csv", row.names = FALSE)
write.csv(summary_table, "original_style_interactions_summary.csv", row.names = FALSE)
write.csv(delta, "original_style_interactions_delta.csv", row.names = FALSE)

plot_metric(results, "y_mean_abs_err", "Absolute Error of mean(y)", "original_style_interactions_abs_err_box.png")
plot_metric(results, "y_mean_bias", "Bias of mean(y)", "original_style_interactions_bias_box.png")
plot_metric(results, "weight_var", "Weight variance", "original_style_interactions_weight_var_box.png")

cat("Completed original-style large-interaction simulation.\n")
cat("Method means:\n")
print(summary_table, row.names = FALSE)
cat("\nPairwise deltas (method_b - method_a):\n")
print(delta, row.names = FALSE)
