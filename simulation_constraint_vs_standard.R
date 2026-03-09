set.seed(20260309)

suppressPackageStartupMessages({
  library(survey)
  library(ggplot2)
})

weighted_cov <- function(x, y, w) {
  w <- w / sum(w)
  mx <- sum(w * x)
  my <- sum(w * y)
  sum(w * (x - mx) * (y - my))
}

weighted_cov_matrix <- function(mat, w) {
  w <- w / sum(w)
  mu <- colSums(mat * w)
  centered <- sweep(mat, 2, mu, "-")
  zw <- centered * sqrt(w)
  crossprod(zw)
}

simulate_population <- function(N_pop = 90000, p_vars = 10, target_rr = 0.35) {
  z <- rnorm(N_pop)

  X <- matrix(0L, nrow = N_pop, ncol = p_vars)
  for (j in seq_len(p_vars)) {
    strength <- runif(1, 0.35, 0.85)
    shift <- runif(1, -0.5, 0.5)
    eta <- shift + strength * z + rnorm(N_pop, sd = 0.8)
    X[, j] <- rbinom(N_pop, 1, plogis(eta))
  }
  colnames(X) <- paste0("X", seq_len(p_vars))

  beta <- seq(0.75, 0.20, length.out = p_vars)
  lin <- as.vector(X %*% beta + 0.20 * z)
  intercept <- uniroot(function(a) mean(plogis(a + lin)) - target_rr,
                       interval = c(-8, 3))$root
  rho <- plogis(intercept + lin)

  e1 <- rnorm(N_pop)
  e2 <- rnorm(N_pop)
  g1 <- seq(0.30, 0.05, length.out = p_vars)
  g2 <- seq(-0.20, 0.20, length.out = p_vars)
  y1 <- 0.9 * z + as.vector(X %*% g1) + e1
  y2 <- 0.6 * z + as.vector(X %*% g2) + 0.7 * e1 + e2

  list(X = X, rho = rho, y1 = y1, y2 = y2)
}

run_simulation <- function(
  reps = 250,
  sample_n = 3200,
  K_values = 2:10,
  N_pop = 90000,
  p_vars = 10,
  target_rr = 0.25
) {
  pop <- simulate_population(N_pop = N_pop, p_vars = p_vars, target_rr = target_rr)
  X <- pop$X
  rho <- pop$rho
  y1 <- pop$y1
  y2 <- pop$y2

  true_cov_y <- cov(y1, y2)
  true_means <- colMeans(X)

  true_cov_X <- lapply(K_values, function(k) cov(X[, seq_len(k), drop = FALSE]))
  names(true_cov_X) <- as.character(K_values)

  out <- vector("list", length = reps * length(K_values) * 2)
  row_id <- 1L

  safe_calibrate <- function(...) {
    fit <- tryCatch(
      withCallingHandlers(
        calibrate(...),
        warning = function(w) invokeRestart("muffleWarning")
      ),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      return(NULL)
    }
    fit
  }

  for (r in seq_len(reps)) {
    sample_idx <- sample.int(N_pop, sample_n, replace = FALSE)
    responded <- rbinom(sample_n, 1, rho[sample_idx]) == 1
    rr_obs <- mean(responded)
    resp_idx <- sample_idx[responded]
    n_resp <- length(resp_idx)

    if (n_resp < 300) {
      next
    }

    df <- as.data.frame(X[resp_idx, , drop = FALSE])
    df$y1 <- y1[resp_idx]
    df$y2 <- y2[resp_idx]
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
        force = TRUE,
        maxit = 300,
        epsilon = 1e-8
      )
      con <- safe_calibrate(
        base_design,
        formula = form,
        population = pop_totals,
        calfun = "raking",
        bounds = c(rr_obs, Inf),
        force = TRUE,
        maxit = 300,
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

        cov_y_hat <- weighted_cov(df$y1, df$y2, w)

        out[[row_id]] <- data.frame(
          rep = r,
          K = k,
          method = method,
          rr_obs = rr_obs,
          weight_var = var(w),
          min_weight = min(w),
          max_weight = max(w),
          cov_mae_x = cov_mae_x,
          cov_abs_err_y = abs(cov_y_hat - true_cov_y)
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
  metrics <- c("weight_var", "cov_mae_x", "cov_abs_err_y", "min_weight", "max_weight")
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
    stat_summary(fun = mean, geom = "point", size = 1.8) +
    labs(x = "Number of weighting variables (K)", y = ylab, color = "") +
    theme_minimal(base_size = 12)
  ggsave(file_name, p, width = 8, height = 5, dpi = 160)
}

results <- run_simulation()
summary_table <- summarize_results(results)

write.csv(results, "constraint_simulation_long.csv", row.names = FALSE)
write.csv(summary_table, "constraint_simulation_summary.csv", row.names = FALSE)

plot_metric(results, "weight_var", "Variance of calibrated weights", "constraint_sim_weight_var.png")
plot_metric(results, "cov_mae_x", "MAE of pairwise covariances among weighting vars", "constraint_sim_cov_mae_x.png")
plot_metric(results, "cov_abs_err_y", "Absolute error of Cov(Y1, Y2)", "constraint_sim_cov_err_y.png")

cat("Completed simulation.\n")
cat(sprintf("Rows in long results: %s\n", nrow(results)))
cat("Summary (means) for key metrics:\n")
print(summary_table[, c("method", "K", "weight_var_mean", "cov_mae_x_mean", "cov_abs_err_y_mean")], row.names = FALSE)
