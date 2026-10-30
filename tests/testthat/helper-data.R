# Simulated data used in several tests.
sim_sparse <- function(n = 300, seed = 1) {
  set.seed(seed)
  x1 <- rexp(n); x2 <- runif(n, 0, 5); x3 <- rnorm(n); x4 <- runif(n)
  u <- runif(n)
  y <- 1 + qnorm(u) + (1 + 2 * u) * x1 + x2
  data.frame(y, x1, x2, x3, x4)
}

# KKT residual of a piqr fit at position l (internal scale), recomputed from scratch.
piqr_kkt <- function(fit, l, data, formula, formula.p) {
  ns <- asNamespace("qrcmNP")
  mf <- model.frame(formula, data = data)
  A <- ns$check_in(mf, formula.p, NULL)
  S <- list(B = A$stats.B, X = A$stats.X, y = A$stats.y)
  # internal theta: invert the (linear) back-transformation numerically
  q <- nrow(A$s); k <- ncol(A$s)
  f <- function(v) c(ns$check_out(matrix(v, q, k), S)$theta)
  b0 <- f(rep(0, q * k))
  Tm <- vapply(seq_len(q * k), function(j) { e <- rep(0, q * k); e[j] <- 1; f(e) - b0 }, numeric(q * k))
  th <- matrix(solve(Tm, c(fit$coefficients[[l]]) - b0), q, k)
  th[abs(th) < 1e-10] <- 0
  ps <- ns$pstar(th, A$V, A$internal.bfun, A$type)
  g <- ns$get_ee(A$type)(th, A$V$y, A$V$z, A$V$d, A$V$X, A$V$Xw, A$internal.bfun, ps$y, ps$z, J = FALSE)$g
  pf <- fit$internal$pf
  free <- A$s == 1
  pen <- list(type = fit$penalty, free = free, W = matrix(pf, q, k) * free, gw = pf * sqrt(rowSums(free)))
  ns$kkt_residual(g, th, fit$lambda[l], pen, which(c(free))) / nrow(data)
}
