test_that("Khatri-Rao product matches the column-major vectorisation of theta", {
  kr <- asNamespace("qrcmNP")$kr
  X <- matrix(1:6, 3, 2); B <- matrix(c(1, 2, 3, 10, 20, 30, 5, 5, 5), 3, 3)
  Z <- kr(X, B)
  theta <- matrix(1:6, 2, 3)
  expect_equal(c(Z %*% c(theta)), rowSums((X %*% theta) * B))
})

test_that("analytic gradient of the integrated loss matches numerical derivatives", {
  ns <- asNamespace("qrcmNP")
  d <- sim_sparse(200)
  mf <- model.frame(y ~ x1 + x2, data = d)
  A <- ns$check_in(mf, ~ slp(p, 2), NULL)
  V <- A$V; bf <- A$internal.bfun
  th <- ns$start_theta(V, bf, A$s, "iqr")
  Lf <- function(t) { t <- matrix(t, nrow(th)); ns$loss_internal(t, V, bf, "iqr", ns$pstar(t, V, bf, "iqr")) }
  ps <- ns$pstar(th, V, bf, "iqr")
  g <- ns$iqr_ee(th, V$y, V$z, V$d, V$X, V$Xw, bf, ps$y, ps$z, J = TRUE)
  ng <- vapply(seq_along(th), function(j) {
    e <- rep(0, length(th)); e[j] <- 1e-3
    (Lf(c(th) + e) - Lf(c(th) - e)) / 2e-3
  }, numeric(1))
  expect_equal(g$g, ng, tolerance = 0.02 * max(abs(ng)), ignore_attr = TRUE)
  expect_equal(g$J, t(g$J))
})

test_that("check_in/check_out round trip and errors", {
  ns <- asNamespace("qrcmNP")
  d <- sim_sparse(100)
  mf <- model.frame(y ~ x1 + x2, data = d)
  expect_error(ns$check_in(mf, ~ slp(p, 2), matrix(1, 2, 2)), "wrong size")
  expect_error(ns$check_in(mf, ~ x1, NULL), "depend on 'p'")
  A <- ns$check_in(mf, ~ slp(p, 2), NULL)
  expect_equal(dim(A$s), c(3, 3))
  expect_true(all(A$V$y >= 0 & A$V$y <= 10))
})
