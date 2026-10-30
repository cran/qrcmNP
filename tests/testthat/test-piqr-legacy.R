# Regression tests against qrcmNP 0.2.2 (CRAN), see data-raw/make_reference_0.2.2.R

test_that("method = 'legacy' reproduces the path of qrcmNP 0.2.2", {
  ref <- readRDS(test_path("fixtures", "ref_0.2.2.rds"))
  r1 <- ref$piqr1
  f1 <- piqr(y ~ x1 + x2, formula.p = ~ I(qnorm(p)) + slp(p, 2), s = r1$s, nlambda = 50,
             data = r1$data, method = "legacy")
  expect_equal(f1$lambda, r1$lambda, tolerance = 1e-3)
  expect_equal(unname(f1$df), unname(r1$df))
  expect_equal(f1$coefficients, r1$coefficients, tolerance = 0.01, ignore_attr = TRUE)
  # the BIC curve is flat at the end of this path: allow a difference of 2 positions
  expect_lte(abs(gof.piqr(f1, method = "BIC", plot = FALSE)$posMinLambda - r1$bic), 2)

  r2 <- ref$piqr2
  f2 <- piqr(y ~ X1 + X2 + X3 + X4 + X5, data = r2$data, nlambda = 30, method = "legacy")
  expect_equal(f2$lambda, r2$lambda, tolerance = 1e-3)
  expect_equal(unname(f2$df), unname(r2$df))
  expect_equal(unname(gof.piqr(f2, method = "BIC", plot = FALSE)$posMinLambda), unname(r2$bic))
  # small differences are due to the starting points (ranks/Kaplan-Meier instead of 'pch')
  expect_equal(f2$obj.function, r2$obj.function, tolerance = 5e-3)
})
