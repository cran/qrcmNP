#' Nonlinear Parametric Quantile Models
#'
#' Estimates a parametric model for the whole conditional quantile function,
#' possibly nonlinear in the parameters, in the covariates, or both, by
#' minimizing the integrated loss function of Frumento and Bottai (2016), as
#' described by Bottai and Cilluffo (2020).
#'
#' @param formula a two-sided formula `y ~ x1 + x2 + ...`; the model matrix
#'   (including the intercept, if any) is passed to `fun2` as `X`.
#' @param data an optional data frame containing the variables in `formula`
#'   (and, possibly, `weights` and `cluster`).
#' @param fun a function `fun(theta, p)` that returns the matrix of coefficient
#'   functions \eqn{\beta(p \mid \theta)}{beta(p | theta)}, with one row for each
#'   element of the vector `p` and one column for each column of the model
#'   matrix. It must be vectorized in `p`.
#' @param fun2 an optional function `fun2(beta, X)` returning the quantile
#'   function \eqn{Q(p_i \mid x_i, \theta)}{Q(p_i | x_i, theta)}, where `beta` is
#'   the matrix returned by `fun` (row \eqn{i} evaluated at \eqn{p_i}) and `X`
#'   the model matrix (row \eqn{i} = \eqn{x_i}). The default (`NULL`) is the model
#'   linear in the covariates, \eqn{Q(p \mid x) = x^\top\beta(p)}{Q(p | x) = x' beta(p)},
#'   i.e. `rowSums(beta * X)`.
#' @param x0 starting values for \eqn{\theta}. If named, the names are used to
#'   label the estimates.
#' @param weights optional non-negative weights: a vector, or a variable in
#'   `data` given without quotes. They are normalized to have mean 1 and are
#'   treated as sampling (inverse-probability) weights in the sandwich formula.
#' @param cluster an optional vector, or a variable in `data` given without
#'   quotes, identifying clusters of correlated observations. It is used to
#'   compute cluster-robust standard errors (sandwich) or a cluster bootstrap.
#' @param control a list of control parameters, see [niqr.control()].
#' @param ... not used. Arguments of the interface of 'qrcmNP' 0.2.x (`X`, `y`)
#'   give an informative error.
#'
#' @details
#' Let \eqn{Q(p \mid x, \theta)}{Q(p | x, theta)} be a parametric model for the
#' conditional quantile function of \eqn{y} given \eqn{x}. The estimator
#' \eqn{\hat\theta} minimizes
#' \deqn{L_n(\theta) = \frac{1}{n}\sum_{i=1}^n w_i \int_0^1
#'   \rho_p\{y_i - Q(p \mid x_i, \theta)\}\,dp,}{%
#'   L(theta) = mean_i w_i int_0^1 rho_p(y_i - Q(p | x_i, theta)) dp,}
#' where \eqn{\rho_p(u) = u\{p - I(u < 0)\}}{rho_p(u) = u (p - I(u < 0))}. The
#' integral is computed with the trapezoidal rule on a grid of `n.p` values of
#' \eqn{p} (denser in the tails), splitting the interval containing
#' \eqn{p^*_i = F(y_i \mid x_i, \theta)}{p*_i = F(y_i | x_i, theta)}, where the
#' integrand is not differentiable. The gradient
#' \eqn{-n^{-1}\sum_i w_i \int_0^1 \{p - I(p \ge p^*_i)\}\,\partial_\theta Q(p \mid x_i, \theta)\,dp}
#' is computed in the same way, with numerical derivatives of `fun` (or of
#' `fun2`) with respect to \eqn{\theta}. The loss is minimized by a damped
#' Newton algorithm that uses the Hessian
#' \deqn{H_n(\theta) = \frac1n\sum_i w_i\,\frac{\partial_\theta Q(p^*_i \mid x_i, \theta)\,
#'   \partial_\theta Q(p^*_i \mid x_i, \theta)^\top}{Q'(p^*_i \mid x_i, \theta)},}{%
#'   H(theta) = mean_i w_i dQ(p*_i) dQ(p*_i)' / Q'(p*_i),}
#' where \eqn{Q'} is the derivative with respect to \eqn{p}, i.e. the inverse of
#' the conditional density. \eqn{H_n} is positive semi-definite and estimates
#' the Hessian of the expected loss at the true parameter.
#'
#' The asymptotic covariance matrix of \eqn{\hat\theta} is estimated by the
#' sandwich formula \eqn{H^{-1}\Omega H^{-1}/n} (Bottai and Cilluffo, 2020,
#' Theorem 3), where \eqn{\Omega} is the covariance of the individual
#' contributions to the gradient (summed within clusters if `cluster` is
#' supplied), or by the (cluster) bootstrap (`niqr.control(se = "boot")`).
#'
#' The model is assumed to be monotonically increasing in \eqn{p} for all
#' observations. A warning is given if quantile crossing is detected.
#' Good starting values `x0` are important: the loss function may have several
#' local minima when \eqn{Q} is not one-to-one in \eqn{\theta}
#' (Bottai and Cilluffo, 2020, Theorem 1).
#'
#' @return An object of class `"niqr"`, a list with items:
#' \item{coefficients}{the estimated \eqn{\theta}.}
#' \item{se}{standard errors.}
#' \item{covar}{the estimated covariance matrix of \eqn{\hat\theta}.}
#' \item{obj.function}{the minimized loss \eqn{L_n(\hat\theta)}.}
#' \item{gradient, hessian, omega}{gradient, Hessian and outer product of the
#'   gradient at \eqn{\hat\theta}.}
#' \item{CDF, PDF}{fitted conditional CDF and density at the observed responses.}
#' \item{converged, n.it}{convergence indicator and number of iterations.}
#' \item{call, terms, X, y, weights, cluster, fun, fun2, control}{the elements
#'   of the model.}
#'
#' @note The matrix interface `niqr(fun, fun2, x0, X, y)` of 'qrcmNP' 0.2.x was
#'   removed in version 1.0.0: use `niqr(y ~ x, data, fun, fun2, x0)`, where the
#'   model matrix of the formula plays the role of the old `X`.
#'
#' @references
#' Bottai, M. and Cilluffo, G. (2020). Nonlinear parametric quantile models.
#' *Statistical Methods in Medical Research*, 29(12), 3757-3769.
#' \doi{10.1177/0962280220941159}
#'
#' Frumento, P. and Bottai, M. (2016). Parametric modeling of quantile
#' regression coefficient functions. *Biometrics*, 72(1), 74-84.
#' \doi{10.1111/biom.12410}
#'
#' @seealso [niqr.control()], [summary.niqr()], [predict.niqr()],
#'   [plot.niqr()], [gof.niqr()], [test.fit.niqr()].
#'
#' @examples
#' ## Model linear in the covariates, nonlinear in the parameters:
#' ## beta0(p) = theta1 + exp(theta2 * p), beta1(p) = theta3 + theta4 * p
#' set.seed(1234)
#' n <- 300
#' x <- runif(n)
#' fun <- function(theta, p) {
#'   cbind(theta[1] + exp(theta[2] * p), theta[3] + theta[4] * p)
#' }
#' u <- runif(n)
#' beta <- fun(c(1, 1, 1, 1), u)
#' y <- beta[, 1] + beta[, 2] * x
#' m <- niqr(y ~ x, fun = fun, x0 = c(a = 0.5, b = 0.5, c = 0.5, d = 0.5))
#' summary(m)
#' predict(m, type = "QF", p = c(0.1, 0.5, 0.9), newdata = data.frame(x = 0.5))
#'
#' ## Model nonlinear in the covariate: Hill dose-response quantile function
#' ## Q(p | dose) = b0(p) + (b1 - b0(p)) / (1 + 10^(b3 - log10(dose)))
#' \donttest{
#' set.seed(1)
#' dose <- rep(10^seq(-2, 2, length.out = 20), each = 15)
#' fun.h <- function(theta, p) {
#'   cbind(theta[1] + exp(theta[2]) * qnorm(p), theta[3], theta[4])
#' }
#' fun2.h <- function(beta, X) {
#'   beta[, 1] + (beta[, 2] - beta[, 1]) / (1 + 10^(beta[, 3] - log10(X[, 2])))
#' }
#' yh <- fun2.h(fun.h(c(1, log(0.5), 10, 0), runif(length(dose))), cbind(1, dose))
#' mh <- niqr(yh ~ dose, fun = fun.h, fun2 = fun2.h, x0 = c(0, 0, 5, 0.5))
#' summary(mh)
#' }
#' @export
niqr <- function(formula, data, fun, fun2 = NULL, x0, weights, cluster,
                 control = niqr.control(), ...) {
  cl <- match.call()
  dots <- names(list(...))
  if (any(c("X", "y") %in% dots) || (!missing(formula) && is.function(formula)))
    stop("the matrix interface niqr(fun, fun2, x0, X, y) was removed in qrcmNP 1.0.0: ",
         "please use niqr(formula, data, fun, fun2, x0), e.g. niqr(y ~ x, data = d, fun = fun, x0 = x0)",
         call. = FALSE)
  if (length(dots)) stop("unused argument(s): ", paste(dots, collapse = ", "), call. = FALSE)
  if (missing(formula) || !inherits(formula, "formula")) stop("'formula' is required", call. = FALSE)
  if (!is.list(control)) stop("'control' must be a list, see niqr.control()", call. = FALSE)
  old <- intersect(names(control), c("display", "n_p", "low_p", "up_p", "alpha", "beta", "maxit_start",
                                     "epsilon", "a1", "h1", "meth", "fun_prime_theta", "fun_prime_beta",
                                     "cluster"))
  if (length(old))
    stop("control parameters of qrcmNP 0.2.x are no longer supported (", paste(old, collapse = ", "),
         "): see ?niqr.control (e.g. 'display' -> 'verbose', 'n_p' -> 'n.p', 'low_p' -> 'eps'; ",
         "'cluster' is now an argument of niqr())", call. = FALSE)
  control <- do.call(niqr.control, control)
  if (missing(fun) || !is.function(fun)) stop("'fun' must be a function of (theta, p)", call. = FALSE)
  if (!is.null(fun2) && !is.function(fun2)) stop("'fun2' must be NULL or a function of (beta, X)", call. = FALSE)
  if (missing(x0)) stop("starting values 'x0' are required", call. = FALSE)

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "weights", "cluster"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  if (is.null(y) || is.matrix(y)) stop("the response must be a numeric vector", call. = FALSE)
  X <- model.matrix(mt, mf)
  xlev <- .getXlevels(mt, mf)
  contr <- attr(X, "contrasts")
  w <- model.weights(mf); if (is.null(w)) w <- rep.int(1, length(y))
  clus <- mf[["(cluster)"]]
  if (length(w) != length(y) || any(w < 0) || anyNA(w)) stop("invalid 'weights'", call. = FALSE)

  fit <- niqr_fit(x0, y, X, w, fun, fun2, control, clus)
  fit$call <- cl
  fit$terms <- mt
  fit$xlevels <- xlev
  fit$contrasts <- contr
  fit
}

#' Control Parameters for niqr
#'
#' @param tol convergence tolerance on the relative change of the parameters
#'   (default `1e-6`).
#' @param maxit maximum number of iterations (default 100).
#' @param n.p number of grid points used to compute the integrals over \eqn{p}
#'   (default 199).
#' @param eps the grid covers \eqn{[\epsilon, 1 - \epsilon]}{[eps, 1 - eps]}
#'   (default `1e-4`).
#' @param h relative step for numerical derivatives (default
#'   `.Machine$double.eps^(1/3)`).
#' @param se how to estimate the covariance matrix: `"sandwich"` (default),
#'   `"boot"` (nonparametric bootstrap, or cluster bootstrap if `cluster` is
#'   supplied), or `"none"`.
#' @param R number of bootstrap replications (default 200).
#' @param verbose logical: print the progress of the algorithm?
#' @return A list of control parameters.
#' @seealso [niqr()]
#' @export
niqr.control <- function(tol = 1e-6, maxit = 100, n.p = 199, eps = 1e-4,
                         h = .Machine$double.eps^(1 / 3), se = c("sandwich", "boot", "none"),
                         R = 200, verbose = FALSE) {
  se <- match.arg(se)
  if (!is.numeric(tol) || tol <= 0) stop("'tol' must be positive", call. = FALSE)
  if (!is.numeric(n.p) || n.p < 10) stop("'n.p' must be at least 10", call. = FALSE)
  if (!is.numeric(eps) || eps <= 0 || eps >= 0.1) stop("'eps' must be in (0, 0.1)", call. = FALSE)
  list(tol = tol, maxit = as.integer(maxit), n.p = as.integer(n.p), eps = eps, h = h,
       se = se, R = as.integer(R), verbose = isTRUE(verbose))
}

#' Fit a niqr model (internal).
#' @noRd
niqr_fit <- function(x0, y, X, w, fun, fun2, control, cluster = NULL, boot = FALSE) {
  theta <- x0
  nm <- names(x0); if (is.null(nm)) nm <- paste0("theta", seq_along(x0) - 1)
  theta <- as.numeric(theta)
  w <- w / mean(w)
  M <- list(X = X, y = y, weights = w, fun = fun, fun2 = fun2, h = control$h,
            grid = niqr_grid(control$n.p, control$eps), n.bisec = 30)
  b0 <- fun(theta, c(0.25, 0.5))
  if (!is.matrix(b0) && length(b0) == 2) b0 <- matrix(b0, nrow = 2)
  if (!is.matrix(b0) || nrow(b0) != 2) stop("'fun(theta, p)' must return a matrix with one row for each value of 'p'", call. = FALSE)
  if (is.null(fun2) && ncol(b0) != ncol(X))
    stop(sprintf("'fun(theta, p)' returns %d columns, but the model matrix has %d columns", ncol(b0), ncol(X)), call. = FALSE)

  ev <- niqr_eval(theta, M)
  if (!is.finite(ev$loss)) stop("the loss function is not finite at 'x0'", call. = FALSE)
  # Levenberg-Marquardt damped Newton: (H + mu * diag(H)) step = -g, with mu
  # increased when the full step does not decrease the loss enough and decreased
  # after a successful step. This prevents huge steps when H is ill-conditioned
  # (e.g., parameters that are weakly identified far from the solution).
  conv <- FALSE
  it <- 0
  mu <- 1e-3
  armijo <- function(th.new, step, g, t) {
    Lnew <- tryCatch(sum(w * niqr_loss_i(th.new, M)) / sum(w), error = function(e) NA)
    list(ok = is.finite(Lnew) && Lnew <= ev$loss + 1e-4 * t * sum(step * g), L = Lnew)
  }
  for (it in seq_len(control$maxit)) {
    # stop when the predicted decrease of a Newton step is negligible (at this
    # point the numerical error of the integrals dominates)
    if (newton_decrement(ev) < 1e-9 * (1 + abs(ev$loss))) { conv <- TRUE; it <- it - 1; break }
    H <- ev$H
    g <- ev$grad
    dH <- pmax(diag(H), 1e-12 * max(1, max(abs(diag(H)))))
    accepted <- FALSE
    for (lm in 1:12) {
      step <- tryCatch(-solve(H + diag(mu * dH, length(g)), g), error = function(e) NULL)
      if (!is.null(step) && all(is.finite(step)) && sum(step * g) < 0) {
        a <- armijo(theta + step, step, g, 1)
        if (a$ok) { accepted <- TRUE; tstep <- 1; break }
      }
      mu <- mu * 10
    }
    fallback <- !accepted
    if (!accepted) {
      # last resort: steepest descent with backtracking
      mu <- 1e-3
      step <- -g / max(dH)
      tstep <- 1
      for (ls in 1:50) {
        a <- armijo(theta + tstep * step, step, g, tstep)
        if (a$ok) { accepted <- TRUE; break }
        tstep <- tstep / 2
      }
    }
    if (!accepted) { conv <- newton_decrement(ev) < 1e-6 * (1 + abs(ev$loss)); break }
    damped <- fallback || mu > 1
    mu <- max(mu / 10, 1e-10)
    change <- max(abs(tstep * step) / pmax(1, abs(theta)))
    theta <- theta + tstep * step
    ev <- niqr_eval(theta, M)
    if (control$verbose) message(sprintf("iter %3d   loss = %.10f   max|grad| = %.3e", it, ev$loss, max(abs(ev$grad))))
    # small steps are a sign of convergence only if they are (nearly) Newton steps
    if (change < control$tol && !damped) { conv <- TRUE; break }
  }
  if (conv) conv <- newton_decrement(ev) < 1e-6 * (1 + abs(ev$loss))
  if (!conv && !boot) warning("the algorithm did not converge", call. = FALSE)

  n <- length(y); d <- length(theta)
  Si <- ev$S * w
  omega <- if (is.null(cluster)) crossprod(Si) / n else {
    Sc <- rowsum(Si, cluster)
    G <- nrow(Sc)
    if (G < 2) stop("at least two clusters are needed", call. = FALSE)
    crossprod(Sc) / n * G / (G - 1)
  }
  covar <- matrix(NA_real_, d, d)
  if (!boot && control$se != "none") {
    Hi <- tryCatch(solve(ev$H), error = function(e) NULL)
    if (is.null(Hi)) warning("the Hessian is singular: standard errors are not available", call. = FALSE)
    else covar <- Hi %*% omega %*% Hi / n
    if (control$se == "boot") covar <- niqr_boot(theta, y, X, w, fun, fun2, control, cluster)
  }
  dimnames(covar) <- list(nm, nm)
  names(theta) <- nm
  se <- sqrt(pmax(diag(covar), 0))
  monotone <- if (!boot) niqr_check_monotone(theta, M) else TRUE
  if (!monotone) warning("quantile crossing detected (the fitted quantile function is decreasing ",
                         "in p for some observations)", call. = FALSE)
  pdf <- ev$pdf
  cdf <- c(ev$ps)

  out <- list(coefficients = theta, se = se, covar = covar, obj.function = ev$loss,
              gradient = ev$grad, hessian = ev$H, omega = omega, CDF = cdf, PDF = pdf,
              converged = conv, n.it = it, X = X, y = y, weights = w, cluster = cluster,
              fun = fun, fun2 = fun2, control = control)
  class(out) <- "niqr"
  out
}

#' Newton decrement g' H^{-1} g (half the predicted decrease of a Newton step).
#' @noRd
newton_decrement <- function(ev) {
  d <- tryCatch(sum(ev$grad * solve(ev$H + diag(1e-10 * max(1, mean(abs(diag(ev$H)))), length(ev$grad)),
                                    ev$grad)), error = function(e) Inf)
  if (!is.finite(d)) sum(ev$grad^2) else abs(d)
}

#' Check that the fitted quantile function is increasing in p on the grid.
#' @noRd
niqr_check_monotone <- function(theta, M) {
  Qg <- niqr_Qgrid(theta, M)
  all(Qg[, -1, drop = FALSE] - Qg[, -ncol(Qg), drop = FALSE] >= -1e-10 * (1 + abs(Qg[, -1, drop = FALSE])))
}

#' Bootstrap covariance matrix (cluster bootstrap if cluster is not NULL).
#' @noRd
niqr_boot <- function(theta, y, X, w, fun, fun2, control, cluster) {
  R <- control$R
  ids <- if (is.null(cluster)) seq_along(y) else cluster
  uid <- unique(ids)
  idx <- split(seq_along(y), factor(ids, levels = uid))
  ctrl <- control; ctrl$se <- "none"; ctrl$verbose <- FALSE
  TH <- matrix(NA_real_, R, length(theta))
  for (r in seq_len(R)) {
    s <- unlist(idx[sample.int(length(uid), length(uid), replace = TRUE)], use.names = FALSE)
    f <- tryCatch(suppressWarnings(niqr_fit(theta, y[s], X[s, , drop = FALSE], w[s], fun, fun2, ctrl,
                                            boot = TRUE)), error = function(e) NULL)
    if (!is.null(f) && f$converged) TH[r, ] <- f$coefficients
  }
  ok <- stats::complete.cases(TH)
  if (sum(ok) < 10) stop("bootstrap failed: too few successful replications", call. = FALSE)
  if (sum(!ok) > 0) warning(sprintf("%d bootstrap replications failed or did not converge (excluded)", sum(!ok)),
                            call. = FALSE)
  stats::cov(TH[ok, , drop = FALSE])
}
