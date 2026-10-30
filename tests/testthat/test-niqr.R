fun <- function(theta, p) cbind(theta[1] + exp(theta[2] * p), theta[3] + theta[4] * p)

test_that("niqr with a model linear in theta reproduces qrcm::iqr", {
  set.seed(3)
  n <- 300; x <- runif(n); u <- runif(n)
  y <- 1 + qnorm(u) + (1 + u) * x
  d <- data.frame(y, x)
  fl <- function(theta, p) cbind(theta[1] + theta[2] * qnorm(p), theta[3] + theta[4] * p)
  m <- niqr(y ~ x, data = d, fun = fl, x0 = c(0, 1, 0, 0))
  mi <- qrcm::iqr(y ~ x, formula.p = ~ I(qnorm(p)) + p, data = d, s = matrix(c(1, 1, 1, 0, 0, 1), 2, 3))
  expect_equal(unname(coef(m)), unname(c(mi$coefficients[1, 1:2], mi$coefficients[2, c(1, 3)])),
               tolerance = 2e-3)
  expect_equal(unname(m$se), unname(sqrt(diag(mi$covar))[c(1, 3, 2, 6)]), tolerance = 0.03)
  expect_equal(m$obj.function * n, c(mi$obj.function), tolerance = 1e-3)
})

test_that("niqr: nonlinear in theta (example of qrcmNP 0.2.2)", {
  ref <- readRDS(test_path("fixtures", "ref_0.2.2.rds"))$niqr1
  m <- niqr(y ~ x, data = ref$data, fun = fun, x0 = rep(0.5, 4))
  expect_true(m$converged)
  expect_equal(unname(coef(m)), ref$coef, tolerance = 0.01)
  expect_lte(m$obj.function, ref$fx + 1e-6)
  expect_true(all(is.finite(m$se)))
})

test_that("the analytic gradient matches numerical derivatives (nonlinear fun2)", {
  ns <- asNamespace("qrcmNP")
  set.seed(1)
  dose <- rep(10^seq(-2, 2, length.out = 10), each = 10)
  f1 <- function(theta, p) cbind(theta[1] + exp(theta[2]) * qnorm(p), theta[3], theta[4])
  f2 <- function(beta, X) beta[, 1] + (beta[, 2] - beta[, 1]) / (1 + 10^(beta[, 3] - log10(X[, 2])))
  y <- f2(f1(c(1, log(0.5), 10, 0), runif(100)), cbind(1, dose))
  M <- list(X = cbind(1, dose), y = y, weights = rep(1, 100), fun = f1, fun2 = f2, h = 6e-6,
            grid = ns$niqr_grid(199, 1e-4), n.bisec = 30)
  th <- c(0.5, 0, 8, 0.3)
  ev <- ns$niqr_eval(th, M)
  Lf <- function(t) mean(ns$niqr_loss_i(t, M))
  ng <- vapply(1:4, function(k) { e <- rep(0, 4); e[k] <- 1e-4; (Lf(th + e) - Lf(th - e)) / 2e-4 }, 1)
  expect_equal(ev$grad, ng, tolerance = 1e-4)
  m <- niqr(y ~ dose, fun = f1, fun2 = f2, x0 = c(0, 0, 5, 0.5))
  expect_equal(unname(coef(m))[3], 10, tolerance = 0.01)
})

test_that("the interfaces removed in 1.0.0 give informative errors", {
  ref <- readRDS(test_path("fixtures", "ref_0.2.2.rds"))$niqr1
  X <- cbind(1, ref$data$x); y <- ref$data$y
  expect_error(niqr(fun = fun, x0 = rep(0.5, 4), X = X, y = y), "removed in qrcmNP 1.0.0")
  expect_error(niqr(fun, NULL, rep(0.5, 4), X, y), "removed in qrcmNP 1.0.0")
  expect_error(niqr(y ~ x, data = ref$data, fun = fun, x0 = rep(0.5, 4),
                    control = list(maxit_start = 5, display = FALSE)), "no longer supported")
  m <- niqr(y ~ x, data = ref$data, fun = fun, x0 = rep(0.5, 4))
  expect_null(m[["x"]]); expect_null(m[["fx"]])
  expect_error(testfit.niqr(m, R = 2), "removed in qrcmNP 1.0.0")
  expect_error(predict(m, type = "QF", p = 0.5, newdata = cbind(x = 0.5)), "data frame")
})

test_that("cluster-robust and bootstrap standard errors", {
  set.seed(2)
  n <- 200; x <- runif(n); id <- rep(1:50, each = 4)
  b <- fun(c(1, 1, 1, 1), runif(n))
  d <- data.frame(y = b[, 1] + b[, 2] * x, x, id)
  m1 <- niqr(y ~ x, data = d, fun = fun, x0 = rep(0.5, 4), cluster = id)
  expect_true(all(m1$se > 0))
  m2 <- niqr(y ~ x, data = d, fun = fun, x0 = rep(0.5, 4), cluster = id,
             control = niqr.control(se = "boot", R = 20))
  expect_true(all(m2$se > 0))
  expect_equal(dim(vcov(m2)), c(4, 4))
})

test_that("niqr methods", {
  ref <- readRDS(test_path("fixtures", "ref_0.2.2.rds"))$niqr1
  m <- niqr(y ~ x, data = ref$data, fun = fun, x0 = c(a = 0.5, b = 0.5, c = 0.5, d = 0.5))
  expect_named(coef(m), c("a", "b", "c", "d"))
  expect_output(print(m), "Coefficients")
  expect_output(print(summary(m)), "AIC")
  expect_output(print(summary(m, p = c(0.1, 0.9))), "p = 0.1")
  b <- predict(m, type = "beta", p = c(0.2, 0.5))
  expect_length(b, 2)
  expect_true(all(c("se", "low", "up") %in% names(b[[1]])))
  qf <- predict(m, type = "QF", p = c(0.1, 0.9), newdata = data.frame(x = c(0, 1)))
  expect_equal(dim(qf$fit), c(2, 2))
  cdf <- predict(m, type = "CDF", newdata = ref$data[1:5, ])
  expect_equal(cdf, m$CDF[1:5], tolerance = 1e-6)
  expect_length(predict(m, type = "sim"), nrow(ref$data))
  g <- gof.niqr(m)
  expect_equal(g$AIC, m$obj.function * nrow(ref$data) + 8)
  tf <- test.fit(m, R = 5)
  expect_equal(dim(tf), c(2, 2))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_silent(plot(m, which = 1:3, ask = FALSE))
})

test_that("niqr input checks", {
  d <- data.frame(y = rnorm(20), x = runif(20))
  expect_error(niqr(y ~ x, data = d, fun = fun), "x0")
  expect_error(niqr(y ~ x, data = d, x0 = 1), "'fun'")
  expect_error(niqr(y ~ x, data = d, fun = function(theta, p) cbind(theta[1] + p), x0 = 1),
               "columns")
})
