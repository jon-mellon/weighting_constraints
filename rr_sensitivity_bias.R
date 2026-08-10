set.seed(20260310)

suppressPackageStartupMessages({
  library(survey)
})

weighted_cov_matrix <- function(mat, w) {
  w <- w / sum(w)
  mu <- colSums(mat * w)
  centered <- sweep(mat, 2, mu, "-")
  zw <- centered * sqrt(w)
  crossprod(zw)
}

simulate_population <- function(N_pop = 90000, p_vars = 25, target_rr = 0.30) {
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

  rr_coefs <- c(
    seq(0.9, 0.35, length.out = 8),
    0.30, 0.25, -0.20, 0.18,
    rep(0.08, p_vars - 12)
  )
  lin_rr <- as.vector(X %*% rr_coefs + 0.30 * z_rr)
  rr_intercept <- uniroot(
    function(a) mean(plogis(a + lin_rr)) - target_rr,
    interval = c(-10, 4)
  )$root
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
    byrow = TRUE,
    ncol = 2
  )
  interaction_coefs <- c(
    1.2, -1.0, 0.9, -0.8, 0.7, -0.6, 0.55,
    -0.45, 0.4, -0.35, 0.3, 0.25, -0.2, 0.18
  )
  interaction_signal <- rowSums(
    vapply(seq_len(nrow(interaction_pairs)), function(ii) {
      p1 <- interaction_pairs[ii, 1]
      p2 <- interaction_pairs[ii, 2]
      interaction_coefs[ii] * (X[, p1] * X[, p2])
    }, numeric(N_pop))
  )

  y <- 0.9 * z_rr + 0.35 * z_aux + as.vector(X %*% y_coefs) +
    interaction_signal + rnorm(N_pop, sd = 1.1)

  list(X = X, rho = rho, y = y, interaction_pairs = interaction_pairs)
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

run_one_rr <- function(
  target_rr,
  reps = 60,
  sample_n = 7000,
  K_values = c(8, 12, 16, 20, 25),
  N_pop = 90000,
  p_vars = 25,
  eval_n = 10000
) {
  pop <- simulate_population(N_pop = N_pop, p_vars = p_vars, target_rr = target_rr)
  X <- pop$X
  rho <- pop$rho
  y <- pop$y
  interaction_pairs <- pop$interaction_pairs

  true_means <- colMeans(X)
  true_y_mean <- mean(y)
  true_cov_X <- lapply(K_values, function(k) cov(X[, seq_len(k), drop = FALSE]))
  names(true_cov_X) <- as.character(K_values)

  eval_idx <- sample.int(N_pop, eval_n, replace = FALSE)
  eval_df <- as.data.frame(X[eval_idx, , drop = FALSE])
  eval_y <- y[eval_idx]

  rows <- vector("list", length = reps * length(K_values) * 2)
  rid <- 1L

  for (r in seq_len(reps)) {
    s_idx <- sample.int(N_pop, sample_n, replace = FALSE)
    resp <- rbinom(sample_n, 1, rho[s_idx]) == 1
    rr_obs <- mean(resp)
    ridx <- s_idx[resp]
    n_resp <- length(ridx)
    if (n_resp < 400) {
      next
    }

    df <- as.data.frame(X[ridx, , drop = FALSE])
    df$y <- y[ridx]
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
        cov_err_vec <- cov_hat[off_diag] - truth[off_diag]
        cov_mae_x <- mean(abs(cov_err_vec))
        cov_bias_x <- mean(cov_err_vec)

        w_norm <- w / sum(w)
        y_mean_hat <- sum(w_norm * df$y)
        y_mean_err <- y_mean_hat - true_y_mean

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
          lm(model_formula, data = df[, c("y", vars), drop = FALSE], weights = w),
          error = function(e) NULL
        )
        y_pred_mae <- NA_real_
        y_pred_bias <- NA_real_
        if (!is.null(fit)) {
          pred <- tryCatch(
            predict(fit, newdata = eval_df[, vars, drop = FALSE]),
            error = function(e) rep(NA_real_, length(eval_y))
          )
          if (!all(is.na(pred))) {
            pe <- pred - eval_y
            y_pred_mae <- mean(abs(pe))
            y_pred_bias <- mean(pe)
          }
        }

        rows[[rid]] <- data.frame(
          target_rr = target_rr,
          rep = r,
          K = k,
          method = method,
          rr_obs = rr_obs,
          weight_var = var(w),
          min_weight = min(w),
          max_weight = max(w),
          cov_mae_x = cov_mae_x,
          cov_bias_x = cov_bias_x,
          y_mean_err = y_mean_err,
          y_mean_abs_err = abs(y_mean_err),
          y_pred_bias = y_pred_bias,
          y_pred_mae = y_pred_mae
        )
        rid <- rid + 1L
      }
    }
  }

  d <- do.call(rbind, rows)
  d <- d[complete.cases(d), ]
  d
}

summarize_bias_mae <- function(d) {
  by <- list(target_rr = d$target_rr, method = d$method, K = d$K)
  agg <- aggregate(
    cbind(
      rr_obs, weight_var, min_weight, max_weight,
      cov_mae_x, cov_bias_x, y_mean_err, y_mean_abs_err, y_pred_bias, y_pred_mae
    ) ~ target_rr + method + K,
    data = d,
    FUN = mean
  )

  # Add variance/RMSE style summaries for y mean and prediction bias
  y_mean_var <- aggregate(y_mean_err ~ target_rr + method + K, data = d, FUN = var)
  names(y_mean_var)[4] <- "y_mean_err_var"
  y_mean_rmse <- aggregate(I(y_mean_err^2) ~ target_rr + method + K, data = d, FUN = function(x) sqrt(mean(x)))
  names(y_mean_rmse)[4] <- "y_mean_rmse"

  y_pred_var <- aggregate(y_pred_bias ~ target_rr + method + K, data = d, FUN = var)
  names(y_pred_var)[4] <- "y_pred_bias_var"
  y_pred_rmse <- aggregate(I(y_pred_bias^2) ~ target_rr + method + K, data = d, FUN = function(x) sqrt(mean(x)))
  names(y_pred_rmse)[4] <- "y_pred_bias_rmse"

  out <- merge(agg, y_mean_var, by = c("target_rr", "method", "K"))
  out <- merge(out, y_mean_rmse, by = c("target_rr", "method", "K"))
  out <- merge(out, y_pred_var, by = c("target_rr", "method", "K"))
  out <- merge(out, y_pred_rmse, by = c("target_rr", "method", "K"))
  out[order(out$target_rr, out$K, out$method), ]
}

make_delta <- function(summary_table) {
  out <- list()
  oi <- 1L
  for (rr in sort(unique(summary_table$target_rr))) {
    for (k in sort(unique(summary_table$K))) {
      s <- subset(summary_table, target_rr == rr & K == k & method == "standard")
      c <- subset(summary_table, target_rr == rr & K == k & method == "constrained")
      if (nrow(s) == 0 || nrow(c) == 0) {
        next
      }
      out[[oi]] <- data.frame(
        target_rr = rr,
        K = k,
        rr_obs = s$rr_obs,
        weight_var_pct = 100 * (c$weight_var - s$weight_var) / s$weight_var,
        cov_mae_x_pct = 100 * (c$cov_mae_x - s$cov_mae_x) / s$cov_mae_x,
        y_mean_abs_err_pct = 100 * (c$y_mean_abs_err - s$y_mean_abs_err) / s$y_mean_abs_err,
        y_pred_mae_pct = 100 * (c$y_pred_mae - s$y_pred_mae) / s$y_pred_mae,
        y_mean_bias_diff = c$y_mean_err - s$y_mean_err,
        y_pred_bias_diff = c$y_pred_bias - s$y_pred_bias,
        y_mean_abs_bias_diff = abs(c$y_mean_err) - abs(s$y_mean_err),
        y_pred_abs_bias_diff = abs(c$y_pred_bias) - abs(s$y_pred_bias),
        min_std = s$min_weight,
        min_con = c$min_weight,
        max_std = s$max_weight,
        max_con = c$max_weight
      )
      oi <- oi + 1L
    }
  }
  do.call(rbind, out)
}

all <- do.call(rbind, lapply(c(0.30, 0.50, 0.65), run_one_rr))
summary_table <- summarize_bias_mae(all)
delta <- make_delta(summary_table)

write.csv(all, "rr_sensitivity_bias_long.csv", row.names = FALSE)
write.csv(summary_table, "rr_sensitivity_bias_summary.csv", row.names = FALSE)
write.csv(delta, "rr_sensitivity_bias_delta.csv", row.names = FALSE)

cat("Done. Counts by target_rr/method/K:\n")
print(with(all, table(target_rr, method, K)))

cat("\nDelta table:\n")
print(delta)

cat("\nAverage delta by target_rr (across available K):\n")
print(aggregate(
  cbind(
    weight_var_pct, cov_mae_x_pct, y_mean_abs_err_pct, y_pred_mae_pct,
    y_mean_bias_diff, y_pred_bias_diff, y_mean_abs_bias_diff, y_pred_abs_bias_diff
  ) ~ target_rr,
  data = delta,
  FUN = mean
))
