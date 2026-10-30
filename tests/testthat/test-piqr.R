d <- sim_sparse(300)
fo <- y ~ x1 + x2 + x3 + x4
fp <- ~ slp(p, 2)

test_that("piqr returns a well-formed path", {
  fit <- piqr(fo, formula.p = fp, data = d, nlambda = 20)
  expect_s3_class(fit, "piqr")
  expect_length(fit$lambda, 20)
  expect_true(all(diff(fit$lambda) < 0))
  expect_length(fit$coefficients, 20)
  expect_equal(dim(fit$coefficients[[1]]), c(5, 3))
  expect_true(all(fit$converged))
  # lambda_max: only the intercept is in the model
  expect_equal(unname(fit$df[1]), 1)
  expect_true(all(fit$coefficients[[1]][-1, ] == 0))
  # the number of parameters is (weakly) increasing along the path at the start
  expect_true(fit$npar[20] >= fit$npar[2])
})

test_that("the end of the path approaches the unpenalized qrcm::iqr fit", {
  fit <- piqr(fo, formula.p = fp, data = d, nlambda = 20, lambda.min.ratio = 1e-6)
  m0 <- qrcm::iqr(fo, formula.p = fp, data = d)
  expect_equal(fit$obj.function[20], c(m0$obj.function), tolerance = 1e-3)
})

test_that("KKT conditions hold along the lasso and group lasso paths", {
  for (pen in c("lasso", "group")) {
    fit <- piqr(fo, formula.p = fp, data = d, nlambda = 10, penalty = pen)
    for (l in c(3, 6, 10)) {
      r <- piqr_kkt(fit, l, d, fo, fp)
      expect_lt(max(abs(r)), 0.02)
    }
  }
})

test_that("group lasso selects whole covariates", {
  fit <- piqr(fo, formula.p = fp, data = d, nlambda = 15, penalty = "group")
  for (th in fit$coefficients) {
    rs <- rowSums(th != 0)
    expect_true(all(rs[-1] %in% c(0, 3)))
  }
  best <- gof.piqr(fit, plot = FALSE)
  sel <- rowSums(fit$seqS[[best$posMinLambda]]) > 0
  expect_equal(unname(sel), c(TRUE, TRUE, TRUE, FALSE, FALSE))
})

test_that("penalty.factor = 0 keeps a covariate unpenalized", {
  fit <- piqr(fo, formula.p = fp, data = d, nlambda = 5, penalty.factor = c(1, 1, 0, 1))
  expect_true(all(vapply(fit$coefficients, function(th) any(th["x3", ] != 0), TRUE)))
  expect_error(piqr(fo, formula.p = fp, data = d, nlambda = 5, penalty.factor = c(1, -1, 0, 1)),
               "non-negative")
  expect_error(piqr(fo, formula.p = fp, data = d, nlambda = 5, penalty.factor = c(1, 1)), "length")
})

test_that("user-supplied lambda and s are respected", {
  s <- matrix(1, 5, 3); s[2:5, 3] <- 0
  fit <- piqr(fo, formula.p = fp, data = d, s = s, lambda = c(10, 1, 0.1))
  expect_equal(fit$lambda, c(10, 1, 0.1))
  expect_true(all(vapply(fit$coefficients, function(th) all(th[2:5, 3] == 0), TRUE)))
  expect_error(piqr(fo, formula.p = fp, data = d, s = matrix(1, 2, 2)), "wrong size")
})

test_that("methods work", {
  fit <- piqr(fo, formula.p = fp, data = d, nlambda = 10)
  expect_error(coef(fit), "pos.lambda")
  expect_equal(coef(fit, pos.lambda = 3), fit$coefficients[[3]])
  expect_equal(coef(fit, lambda = fit$lambda[4]), fit$coefficients[[4]])
  expect_equal(nobs(fit), 300)
  expect_s3_class(summary(fit, pos.lambda = 5), "summary.piqr")
  expect_output(print(summary(fit, pos.lambda = 5)), "Selected value of lambda")
  expect_output(print(summary(fit, pos.lambda = 5, p = c(0.1, 0.5))), "p = 0.1")
  expect_output(print(fit), "Penalty")
  b <- predict(fit, pos.lambda = 10, type = "beta", p = c(0.2, 0.8))
  expect_length(b, 5)
  qf <- predict(fit, pos.lambda = 10, type = "QF", p = c(0.25, 0.5), newdata = d[1:3, ])
  expect_equal(dim(qf), c(3, 2))
  cdf <- predict(fit, pos.lambda = 10, type = "CDF")
  expect_true(all(cdf$CDF >= 0 & cdf$CDF <= 1))
  expect_length(predict(fit, pos.lambda = 10, type = "sim"), 300)
  g <- gof.piqr(fit, plot = FALSE)
  expect_s3_class(g, "gof.piqr")
  expect_output(print(g), "Selected lambda")
  g2 <- gof.piqr(fit, Cn = "log(log(n))", plot = FALSE, df = "covariates")
  expect_true(g2$posMinLambda >= 1)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_silent(plot(fit, xvar = "lambda"))
  expect_silent(plot(fit, xvar = "grad"))
  expect_silent(plot(fit, xvar = "objective"))
  expect_silent(plot(fit, xvar = "beta", pos.lambda = 10, which = 1:6, ask = FALSE))
})

test_that("refit_piqr re-estimates the selected model with qrcm::iqr", {
  fit <- piqr(fo, formula.p = fp, data = d, nlambda = 15)
  m <- refit_piqr(fit)
  expect_s3_class(m, "iqr")
  sel <- rowSums(fit$seqS[[attr(m, "pos.lambda")]]) > 0
  expect_true(all(m$coefficients[!sel, ] == 0))
  expect_true(all(sqrt(diag(m$covar)) >= 0))
})

test_that("removed arguments give informative errors", {
  expect_error(piqr(fo, formula.p = fp, data = d, nlambda = 2, display = FALSE), "removed in qrcmNP 1.0.0")
  expect_error(piqr(fo, formula.p = fp, data = d, nlambda = 2, foo = 1), "unused argument")
  fit <- piqr(fo, formula.p = fp, data = d, nlambda = 3)
  expect_error(gof.piqr(fit, df.new = TRUE, plot = FALSE), "removed in qrcmNP 1.0.0")
})

test_that("censored data are supported", {
  set.seed(2)
  cens <- rexp(nrow(d), 0.05)
  dc <- data.frame(t = pmin(d$y + 5, cens), status = as.numeric(d$y + 5 <= cens), d[, -1])
  fl <- piqr(survival::Surv(t, status) ~ x1 + x2 + x3, formula.p = fp, data = dc, nlambda = 10)
  fg <- piqr(survival::Surv(t, status) ~ x1 + x2 + x3, formula.p = fp, data = dc, nlambda = 10,
             penalty = "group", lambda = fl$lambda)
  expect_true(all(fl$converged))
  expect_true(all(is.finite(fl$obj.function)))
  expect_equal(fl$coefficients[[10]], fg$coefficients[[10]], tolerance = 0.05)
  expect_s3_class(gof.piqr(fl, plot = FALSE), "gof.piqr")
})

test_that("cv_piqr works", {
  skip_on_cran()
  cv <- cv_piqr(fo, formula.p = fp, data = d, nfolds = 3, nlambda = 10)
  expect_s3_class(cv, "cv.piqr")
  expect_length(cv$cvm, 10)
  expect_true(cv$pos.1se <= cv$pos.min)
  expect_output(print(cv), "Cross-validated")
})
