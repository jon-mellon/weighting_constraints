set.seed(20260309)
suppressPackageStartupMessages(library(survey))

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
  beta <- seq(0.95, 0.30, length.out = p_vars)
  lin <- as.vector(X %*% beta + 0.30 * z)
  intercept <- uniroot(function(a) mean(plogis(a + lin)) - target_rr,
                       interval = c(-8, 3))$root
  rho <- plogis(intercept + lin)
  list(X = X, rho = rho)
}

pop <- simulate_population()
X <- pop$X
rho <- pop$rho
means <- colMeans(X)
print(means)

check_cal <- function(des, form, pop_totals, rr) {
  msg_std <- NULL
  msg_con <- NULL
  std <- withCallingHandlers(
    tryCatch(calibrate(des, formula = form, population = pop_totals, calfun = "raking", maxit = 300, epsilon = 1e-8),
             error = function(e) e),
    warning = function(w) {
      msg_std <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )

  con <- withCallingHandlers(
    tryCatch(calibrate(des, formula = form, population = pop_totals, calfun = "raking", bounds = c(rr, Inf), maxit = 300, epsilon = 1e-8),
             error = function(e) e),
    warning = function(w) {
      msg_con <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )

  list(std_msg = msg_std, con_msg = msg_con,
       std_err = inherits(std, "error"), con_err = inherits(con, "error"))
}

sample_n <- 3200
for (r in 1:20) {
  s <- sample.int(nrow(X), sample_n, replace = FALSE)
  resp <- rbinom(sample_n,1,rho[s])==1
  rr <- mean(resp)
  ridx <- s[resp]
  df <- as.data.frame(X[ridx,,drop=FALSE])
  des <- svydesign(ids=~1,data=df,weights=~1)

  for (k in 2:5) {
    vars <- paste0("X",1:k)
    form <- as.formula(paste("~",paste(vars,collapse=" + ")))
    pop_totals <- c("(Intercept)"=nrow(df), setNames(nrow(df)*means[1:k], vars))
    chk <- check_cal(des, form, pop_totals, rr)
    if (!is.null(chk$std_msg) || !is.null(chk$con_msg) || chk$std_err || chk$con_err) {
      cat("rep",r,"k",k,"rr",round(rr,3),"std_msg",chk$std_msg,"con_msg",chk$con_msg,
          "std_err",chk$std_err,"con_err",chk$con_err,"\n")
    }
  }
}
