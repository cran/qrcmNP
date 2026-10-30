# Edge cases (see the code audit in NOTES/qrcmNP)

test_that("piqr without unpenalized coefficients (no intercept)", {
  d <- sim_sparse(200)
  fit <- piqr(y ~ 0 + x1 + x2 + x3, formula.p = ~ slp(p, 2), data = d, nlambda = 8)
  expect_true(all(fit$converged))
  expect_true(all(vapply(fit$coefficients, function(th) all(is.finite(th)), TRUE)))
  expect_true(fit$npar[8] > 0)
})

test_that("piqr with factors, weights and missing values", {
  d <- sim_sparse(240)
  d$f <- factor(rep(c("a", "b", "c"), 80))
  d$y <- d$y + (d$f == "b")
  d$w <- rep(c(1, 2), 120)
  d$x3[5] <- NA
  fit <- piqr(y ~ x1 + x2 + x3 + f, formula.p = ~ slp(p, 2), data = d, weights = w, nlambda = 8)
  expect_equal(nobs(fit), 239)
  expect_equal(nrow(fit$coefficients[[1]]), 6)
  qf <- predict(fit, pos.lambda = 8, type = "QF", p = 0.5, newdata = d[1:6, ])
  expect_equal(nrow(qf), 6)
  expect_true(is.na(qf[5, 1]))
  cdf <- predict(fit, pos.lambda = 8, type = "CDF", newdata = d[1:6, ])
  expect_equal(nrow(cdf), 6)
  cv <- cv_piqr(y ~ x1 + x2 + x3 + f, formula.p = ~ slp(p, 2), data = d, weights = w,
                nlambda = 8, nfolds = 3)
  expect_length(cv$cvm, 8)
  expect_true(all(is.finite(cv$cvm)))
  expect_false(any(grepl("data.frame", deparse(cv$fit$call))))
})

test_that("a single (user-supplied) lambda", {
  d <- sim_sparse(200)
  fit <- piqr(y ~ x1 + x2 + x3, formula.p = ~ slp(p, 2), data = d, lambda = 1)
  expect_length(fit$lambda, 1)
  expect_equal(coef(fit), fit$coefficients[[1]])
  expect_s3_class(summary(fit), "summary.piqr")
  fit1 <- piqr(y ~ x1 + x2 + x3, formula.p = ~ slp(p, 2), data = d, nlambda = 1)
  expect_equal(unname(fit1$df), 1)
})

test_that("non-slp basis: selected coefficients are counted on the internal scale", {
  d <- sim_sparse(200)
  fit <- piqr(y ~ x1 + x2 + x3, formula.p = ~ I(qnorm(p)) + p, data = d, nlambda = 8)
  expect_equal(unname(fit$npar), unname(vapply(fit$seqS, sum, 1)))
  expect_equal(unname(fit$npar[1]), 3)
})

test_that("niqr with a vector-valued fun, q = 1", {
  set.seed(1)
  y <- 2 + 0.5 * qnorm(runif(200))
  m <- niqr(y ~ 1, fun = function(theta, p) theta[1] + exp(theta[2]) * qnorm(p), x0 = c(0, 0))
  expect_true(m$converged)
  expect_equal(unname(coef(m)), c(2, log(0.5)), tolerance = 0.15)
  s <- summary(m, p = 0.5)
  expect_s3_class(s, "summary.niqr")
})

test_that("niqr predictions with NA in newdata, and no standard errors", {
  ref <- readRDS(test_path("fixtures", "ref_0.2.2.rds"))$niqr1
  fun <- function(theta, p) cbind(theta[1] + exp(theta[2] * p), theta[3] + theta[4] * p)
  m <- niqr(y ~ x, data = ref$data, fun = fun, x0 = rep(0.5, 4), control = niqr.control(se = "none"))
  nd <- data.frame(x = c(0.1, NA, 0.9))
  qf <- predict(m, type = "QF", p = 0.5, newdata = nd, se = FALSE)
  expect_equal(nrow(qf), 3)
  expect_true(is.na(qf[2, 1]))
  expect_output(print(summary(m, p = 0.5)), "p = 0.5")
  expect_error(niqr(y ~ x, data = ref$data, fun = fun, x0 = rep(0.5, 4), cluster = rep(1, 300)),
               "two clusters")
})
