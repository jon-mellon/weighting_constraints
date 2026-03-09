set.seed(1)
library(survey)

N <- 50000
P <- 8
z <- rnorm(N)
X <- sapply(1:P, function(j) {
  rbinom(N, 1, plogis(-0.1 + runif(1, 0.4, 0.9) * z + rnorm(N, sd = 0.7)))
})
colnames(X) <- paste0("X", 1:P)

beta <- seq(0.9, 0.3, length.out = P)
lin <- as.vector(X %*% beta + 0.2 * z)
int <- uniroot(function(a) mean(plogis(a + lin)) - 0.35, c(-8, 2))$root
rho <- plogis(int + lin)

n <- 3000
s <- sample.int(N, n)
R <- rbinom(n, 1, rho[s])
rr <- mean(R)
ridx <- s[R == 1]

df <- as.data.frame(X[ridx, , drop = FALSE])
k <- 8
vars <- paste0("X", 1:k)
mm <- colMeans(X[, 1:k, drop = FALSE])

des <- svydesign(ids = ~1, data = df, weights = ~1)
form <- as.formula(paste("~", paste(vars, collapse = "+")))
poptot <- c("(Intercept)" = nrow(df), setNames(nrow(df) * mm, vars))

std <- calibrate(design = des, formula = form, population = poptot, calfun = "raking")
con <- calibrate(design = des, formula = form, population = poptot, calfun = "raking", bounds = c(rr, Inf))

w1 <- weights(std)
w2 <- weights(con)
cat("rr", rr, "\n")
cat("std min max var", min(w1), max(w1), var(w1), "\n")
cat("con min max var", min(w2), max(w2), var(w2), "\n")
