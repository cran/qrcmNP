# qrcmNP 1.0.0

## Breaking changes

* `piqr()` now solves by default the lasso problem of Sottile et al. (2020, eq. 3.4)
  (`method = "pnewton"`): results differ from version 0.2.x. The previous algorithm
  is available with `method = "legacy"`.
* `library(qrcmNP)` no longer attaches 'qrcm' and 'survival' (moved to `Imports`);
  `slp()`, `plf()` and `test.fit()` are re-exported. Use `library(qrcm)`,
  `library(survival)` or `qrcm::iqr()`, `survival::Surv()` where needed.
* `gof.piqr()`: the default `df` counts the nonzero coefficients (`"parameters"`);
  use `df = "covariates"` for the previous definition (default with `method = "legacy"`).
* Default output is silent (`verbose = FALSE`); `ask` defaults to `interactive()`.
* `niqr()`: standard errors are now computed with the analytic Hessian, so they
  differ from those of version 0.2.x.
* Removed (with an informative error when used):

  | Removed in 1.0.0 | Use instead |
  |---|---|
  | `testfit.niqr(obj, R)` | `test.fit(object, R)` |
  | `niqr(fun, fun2, x0, X, y, control)` | `niqr(formula, data, fun, fun2, x0)` |
  | `piqr(..., display = )` | `piqr(..., verbose = )` |
  | `gof.piqr(..., df.new = )` | `gof.piqr(..., df = )` |
  | `niqr()` control names `display`, `n_p`, `low_p`, `up_p`, `maxit_start`, `h1`, `meth`, ... | `niqr.control(verbose, n.p, eps, ...)`; `cluster` is an argument of `niqr()` |
  | components `x`, `fx`, `gx`, `Hx`, `p.star.y`, `code`, `iter` of `niqr` objects | `coefficients`, `obj.function`, `gradient`, `hessian`, `CDF`, `converged`, `n.it` |
  | matrix `newdata` (with a column `"y"`) in `predict.niqr()` | a data frame |


## Penalized models (`piqr()`)

* New solver (`method = "pnewton"`, default): proximal Newton algorithm with exact
  gradient and Hessian of the integrated loss, coordinate descent for the penalized
  quadratic subproblem, warm starts, sequential strong rules and KKT checks. It
  solves exactly the penalized problem of Sottile et al. (2020, eq. 3.4).
* New argument `penalty = c("lasso", "group")`: the group lasso selects covariates
  as a whole (all the coefficients of a covariate are zero or nonzero).
* New argument `penalty.factor` (covariate-specific penalty weights; 0 = unpenalized).
* The algorithm of version 0.2.x is still available with `method = "legacy"`.
* New `refit_piqr()`: unpenalized re-estimation of the selected model with
  `qrcm::iqr()` (post-selection inference, Sottile et al., 2020, Section 3.2).
* New `cv_piqr()`: K-fold cross-validation of the integrated loss.
* `gof.piqr()`: `Cn` can be a number; new argument `df` (`"parameters"`, the default
  for the new solver, or `"covariates"`, the default for `method = "legacy"`); the
  result has a `print()` method and contains the criterion along the path.
* `obj.function` is now computed also for censored and truncated data, so that
  `gof.piqr()` works with `Surv()` responses.
* `coef()`, `summary()` and `predict()` gain the argument `lambda`; they no longer
  fail when the path contains a single value of lambda.
* `plot(xvar = "beta")`: the QQ-plot now uses `ppoints()` (it previously used sorted
  uniform random numbers); `ask` defaults to `interactive()`.
* Messages are issued with `message()` (argument `verbose`).
* Covariance matrices are computed only by `method = "legacy"`; use `refit_piqr()`.

## Nonlinear models (`niqr()`)

* New formula interface: `niqr(formula, data, fun, fun2, x0, weights, cluster, control)`.
* The gradient is now correct also when `fun2` is nonlinear in `beta` or in the
  covariates (e.g., Hill dose-response models); it was correct only for models
  linear in the covariates.
* Integrals over p are computed on a grid that is denser in the tails, splitting the
  cell that contains F(y | x), which makes loss and gradient smooth and consistent.
  The loss no longer requires O(n^2) operations.
* Optimization by Levenberg-Marquardt damped Newton with the analytic Gauss-Newton
  Hessian, instead of BFGS.
* Standard errors: sandwich formula with the analytic Hessian (previously the BFGS
  approximation, or the identity matrix if the algorithm stopped within 10 iterations),
  cluster-robust sandwich, bootstrap and cluster bootstrap (`niqr.control(se = "boot")`).
* New `niqr.control()`.
* New methods `coef()`, `vcov()`, `nobs()`, `terms()`, `model.matrix()`; `predict()`
  gains `type = "PDF"` and accepts data frames as `newdata`; `summary()` reports the
  information criteria of Bottai and Cilluffo (2020).
* New `gof.niqr()` (AIC, BIC and their log versions, for one or more models).
* `test.fit()`: the Kolmogorov-Smirnov statistic is now two-sided.
* New `test.fit()` method for `niqr` objects; failed bootstrap replications no longer
  stop the procedure.

## Package infrastructure

* The package no longer calls non-exported functions of 'qrcm' through
  `getFromNamespace()`: the required internal routines (adapted from 'qrcm', with
  attribution) are now part of 'qrcmNP'. Jacobians are computed with a single
  cross-product of Khatri-Rao matrices.
* 'qrcm' and 'survival' moved from `Depends` to `Imports`; `slp()`, `plf()` and
  `test.fit()` are re-exported from 'qrcm'. Note that `library(qrcmNP)` no longer
  attaches 'qrcm' and 'survival': scripts calling, e.g., `iqr()` or `Surv()` need
  `library(qrcm)` / `library(survival)` (or `qrcm::iqr()`, `survival::Surv()`).
* Documentation written with roxygen2; new vignette, tests (testthat) and NEWS.

# qrcmNP 0.2.2

* Included the Authors@R field in the DESCRIPTION.
* Minor fixes to the linking targets problem.

# qrcmNP 0.2.0

* Renamed internal functions.
* Improved starting points.
* New dependency 'pch' (>= 2.0).
