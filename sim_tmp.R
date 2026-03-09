set.seed(20260309)

weighted_cov <- function(x, y, w) {
  w <- w / sum(w)
  mx <- sum(w * x)
  my <- sum(w * y)
  sum(w * (x - mx) * (y - my))
}

weighted_cov_matrix <- function(mat, w) {
  w <- w / sum(w)
  mu <- colSums(mat * w)
  xc <- sweep(mat, 2, mu, "-")
  zw <- xc * sqrt(w)
  crossprod(zw)
}

rake_binary <- function(df, targets, vars, min_wt = NULL, max_iter = 200, tol = 1e-7) {
  n <- nrow(df)
  w <- rep(1, n)

  for (iter in seq_len(max_iter)) {
    w_prev <- w

    for (v in vars) {
      x <- df[[v]]
      tot <- sum(w)
      p1 <- sum(w[x == 1]) / tot
      p0 <- 1 - p1
      t1 <- targets[[v]]
      t0 <- 1 - t1

      if (p1 > 0) w[x == 1] <- w[x == 1] * (t1 / p1)
      if (p0 > 0) w[x == 0] <- w[x == 0] * (t0 / p0)
    }

    w <- w / mean(w)

    if (!is.null(min_wt) && any(w < min_wt)) {
      low <- w < min_wt
      low_total <- sum(rep(min_wt, sum(low)))
      high_total <- sum(w[!low])
      remaining <- length(w) - low_total

      if (remaining > 0 && high_total > 0) {
        scale_high <- remaining / high_total
        w[low] <- min_wt
        w[!low] <- w[!low] * scale_high
      } else {
        w[] <- min_wt
      }
    }

    errs <- vapply(vars, function(v) {
      x <- df[[v]]
      abs(sum(w[x == 1]) / sum(w) - targets[[v]])
    }, numeric(1))

    if (max(errs) < tol && max(abs(w - w_prev)) < tol) {
      break
    }
  }

  list(weights = w, iterations = iter, max_margin_error = max(errs))
}

# Population setup
N <- 80000
P <- 10
K_values <- 2:10
R <- 200
n_sample <- 3000

z <- rnorm(N)
X <- matrix(0L, nrow = N, ncol = P)
for (j in seq_len(P)) {
  a <- runif(1, 0.35, 0.9)
  b <- runif(1, -0.4, 0.4)
  eta <- b + a * z + rnorm(N, sd = 0.7)
  X[, j] <- rbinom(N, 1, plogis(eta))
}
colnames(X) <- paste0("X", seq_len(P))

beta <- seq(0.9, 0.25, length.out = P)
lin_base <- as.vector(X %*% beta + 0.25 * z)
# Solve intercept for target overall RR ~ 0.35
f <- function(a) mean(plogis(a + lin_base)) - 0.35
int <- uniroot(f, interval = c(-8, 2))$root
rho <- plogis(int + lin_base)

# Analysis vars
u <- rnorm(N)
e1 <- rnorm(N)
e2 <- rnorm(N)
g1 <- seq(0.35, 0.05, length.out = P)
g2 <- seq(-0.25, 0.2, length.out = P)
y1 <- 0.8 * z + as.vector(X %*% g1) + e1
y2 <- 0.55 * z + as.vector(X %*% g2) + 0.7 * e1 + e2 + 0.25 * u

true_cov_y <- cov(y1, y2)
true_cov_X <- lapply(K_values, function(k) cov(X[, seq_len(k), drop = FALSE]))
names(true_cov_X) <- as.character(K_values)
true_targets <- lapply(K_values, function(k) colMeans(X[, seq_len(k), drop = FALSE]))
names(true_targets) <- as.character(K_values)

cat(sprintf("Population RR target achieved: %.3f\n", mean(rho)))

results <- vector("list", length = R * length(K_values) * 2)
idx <- 1L

for (r in seq_len(R)) {
  s_idx <- sample.int(N, n_sample, replace = FALSE)
  resp <- rbinom(n_sample, 1, rho[s_idx]) == 1
  rr_obs <- mean(resp)
  ridx <- s_idx[resp]

  if (length(ridx) < 300) next

  df <- as.data.frame(X[ridx, , drop = FALSE])
  df$y1 <- y1[ridx]
  df$y2 <- y2[ridx]

  for (k in K_values) {
    vars <- paste0("X", seq_len(k))
    targets <- as.list(true_targets[[as.character(k)]])
    names(targets) <- vars

    std <- rake_binary(df = df, targets = targets, vars = vars, min_wt = NULL, max_iter = 250, tol = 1e-6)
    con <- rake_binary(df = df, targets = targets, vars = vars, min_wt = rr_obs, max_iter = 250, tol = 1e-6)

    for (meth in c("standard", "constrained")) {
      out <- if (meth == "standard") std else con
      w <- out$weights
      cov_y_hat <- weighted_cov(df$y1, df$y2, w)
      cov_mat <- weighted_cov_matrix(as.matrix(df[, vars, drop = FALSE]), w)
      truth <- true_cov_X[[as.character(k)]]
      cov_err <- mean(abs(cov_mat[upper.tri(cov_mat)] - truth[upper.tri(truth)]))

      results[[idx]] <- data.frame(
        rep = r,
        K = k,
        method = meth,
        rr_obs = rr_obs,
        weight_var = var(w),
        min_weight = min(w),
        max_weight = max(w),
        margin_error = out$max_margin_error,
        cov_mae_X = cov_err,
        cov_abs_err_y = abs(cov_y_hat - true_cov_y)
      )
      idx <- idx + 1L
    }
  }
}

res <- do.call(rbind, results)
res <- res[complete.cases(res), ]

aggregate_stats <- function(df, metric) {
  do.call(rbind, lapply(split(df, list(df$method, df$K), drop = TRUE), function(s) {
    data.frame(
      method = s$method[1],
      K = s$K[1],
      n = nrow(s),
      mean = mean(s[[metric]]),
      sd = sd(s[[metric]])
    )
  }))
}

wv <- aggregate_stats(res, "weight_var")
cx <- aggregate_stats(res, "cov_mae_X")
cy <- aggregate_stats(res, "cov_abs_err_y")
me <- aggregate_stats(res, "margin_error")

summ <- merge(merge(wv[, c("method", "K", "mean")], cx[, c("method", "K", "mean")], by = c("method", "K"), suffixes = c("_weight_var", "_cov_mae_X")),
              cy[, c("method", "K", "mean")], by = c("method", "K"))
names(summ)[names(summ) == "mean"] <- "mean_cov_abs_err_y"

print(head(res))
cat("\nMean margin error by method/K (should be tiny):\n")
print(me)
cat("\nSummary means:\n")
print(summ[order(summ$K, summ$method), ])

write.csv(res, "simulation_results_long.csv", row.names = FALSE)
write.csv(summ[order(summ$K, summ$method), ], "simulation_results_summary.csv", row.names = FALSE)
