set.seed(20260310)

suppressPackageStartupMessages({
  library(survey)
  library(ggplot2)
})

# This scenario is tuned so that a real lower-bound constraint on weights
# helps recover the population mean of y under omitted interaction structure.
simulate_population <- function(N_pop = 90000) {
  p <- 14
  z1 <- rnorm(N_pop)
  z2 <- rnorm(N_pop)

  X <- matrix(0L, nrow = N_pop, ncol = p)
  for (j in seq_len(p)) {
    eta <- runif(1, 0.4, 1.1) * z1 + runif(1, -0.5, 0.5) * z2 + rnorm(N_pop, sd = 1.0)
    X[, j] <- rbinom(N_pop, 1, plogis(eta + runif(1, -0.4, 0.4)))
  }
  colnames(X) <- paste0("X", seq_len(p))

  # Hidden response class driven by interactions.
  H <- as.integer((X[, 1] * X[, 2] + X[, 3] * X[, 4] + X[, 1] * X[, 5]) >= 2)

  target_rr <- 0.45
  base_rr <- 0.12
  resp_int_strength <- 1.6

  lin <- qlogis(base_rr) +
    0.6 * X[, 1] + 0.5 * X[, 2] - 0.4 * X[, 6] + 0.3 * X[, 7] +
    resp_int_strength * (X[, 1] * X[, 2]) - 0.8 * (X[, 3] * X[, 4]) +
    0.6 * z1
  rho <- plogis(lin)
  rho <- ifelse(H == 1, 1, rho)

  # Rescale non-H response propensities to hit target overall RR.
  for (ii in seq_len(40)) {
    mult <- target_rr / mean(rho)
    rho <- ifelse(H == 1, 1, pmin(pmax(rho * mult, 0.01), 0.98))
  }

  y_bonus <- -2.0
  y_noise <- 1.6
  y_int_strength <- 0.6
  y <- 0.7 * z1 + 0.3 * z2 + 0.6 * X[, 1] + 0.45 * X[, 2] - 0.35 * X[, 6] + 0.25 * X[, 7] +
    y_int_strength * (X[, 1] * X[, 2]) - 0.9 * (X[, 3] * X[, 4]) +
    y_bonus * H + rnorm(N_pop, sd = y_noise)

  list(X = X, y = y, rho = rho, H = H)
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

run_simulation <- function(reps = 300, n_selected = 1800, K = 10) {
  pop <- simulate_population()
  X <- pop$X
  y <- pop$y
  rho <- pop$rho

  true_y_mean <- mean(y)
  true_means <- colMeans(X)

  vars <- paste0("X", seq_len(K))
  form <- as.formula(paste("~", paste(vars, collapse = " + ")))

  out <- vector("list", length = reps * 3)
  row_id <- 1L

  for (r in seq_len(reps)) {
    s_idx <- sample.int(nrow(X), n_selected, replace = FALSE)
    responded <- rbinom(n_selected, 1, rho[s_idx]) == 1
    rr_obs <- mean(responded)
    ridx <- s_idx[responded]

    if (length(ridx) < 300) {
      next
    }

    df <- as.data.frame(X[ridx, , drop = FALSE])
    df$y <- y[ridx]

    des <- svydesign(ids = ~1, data = df, weights = ~1)
    pop_totals <- c("(Intercept)" = nrow(df), setNames(nrow(df) * true_means[seq_len(K)], vars))

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

    oracle_w <- (1 / rho[ridx])
    oracle_w <- oracle_w / mean(oracle_w)

    for (method in c("standard", "constrained", "oracle_ipw")) {
      w <- if (method == "standard") {
        as.numeric(weights(std))
      } else if (method == "constrained") {
        as.numeric(weights(con))
      } else {
        oracle_w
      }
      wn <- w / sum(w)
      y_hat <- sum(wn * df$y)

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
    bias_diff_b_minus_a = mean(paired$y_mean_bias_b - paired$y_mean_bias_a, na.rm = TRUE),
    abs_err_diff_b_minus_a = mean(paired$y_mean_abs_err_b - paired$y_mean_abs_err_a, na.rm = TRUE),
    weight_var_pct_diff_b_minus_a = mean(100 * (paired$weight_var_b - paired$weight_var_a) / paired$weight_var_a, na.rm = TRUE),
    share_abs_err_improved_b_vs_a = mean(paired$y_mean_abs_err_b < paired$y_mean_abs_err_a, na.rm = TRUE)
  )
}

plot_metric <- function(res, metric, ylab, path) {
  p <- ggplot(res, aes(x = method, y = .data[[metric]], fill = method)) +
    geom_boxplot(width = 0.6, outlier.alpha = 0.25) +
    labs(x = "", y = ylab) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
  ggsave(path, p, width = 6.5, height = 4.5, dpi = 160)
}

results <- run_simulation()
summary_table <- summarize_results(results)

delta <- rbind(
  compare_methods(results, "standard", "constrained"),
  compare_methods(results, "standard", "oracle_ipw"),
  compare_methods(results, "constrained", "oracle_ipw")
)

write.csv(results, "constraint_helps_long.csv", row.names = FALSE)
write.csv(summary_table, "constraint_helps_summary.csv", row.names = FALSE)
write.csv(delta, "constraint_helps_delta.csv", row.names = FALSE)

plot_metric(results, "y_mean_abs_err", "Absolute Error of mean(y)", "constraint_helps_abs_err_box.png")
plot_metric(results, "y_mean_bias", "Bias of mean(y)", "constraint_helps_bias_box.png")
plot_metric(results, "weight_var", "Weight variance", "constraint_helps_weight_var_box.png")

cat("Completed constraint-helps simulation.\n")
cat("Method means:\n")
print(summary_table, row.names = FALSE)
cat("\nPairwise deltas (method_b - method_a):\n")
print(delta, row.names = FALSE)
