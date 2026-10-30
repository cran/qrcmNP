#' Penalized Quantile Regression Coefficients Modeling
#'
#' Fits the penalized integrated loss minimization (PILM) estimator of
#' Sottile et al. (2020) along a path of values of the tuning parameter
#' \eqn{\lambda}. Quantile regression coefficients are modelled as parametric
#' functions of the quantile order \eqn{p}, and covariates are selected using
#' the information on all quantiles simultaneously.
#'
#' @param formula a two-sided formula of the form `y ~ x1 + x2 + ...`. The
#'   response can be a numeric vector or a [survival::Surv()] object of type
#'   `"right"` (right censoring) or `"counting"` (left truncation and right
#'   censoring), as in [qrcm::iqr()].
#' @param formula.p a one-sided formula of the form `~ b1(p, ...) + b2(p, ...) + ...`
#'   describing how the quantile regression coefficients depend on `p`. The
#'   default is a shifted Legendre polynomial of degree 3, see [qrcm::slp()].
#' @param weights an optional vector of non-negative weights.
#' @param data an optional data frame, list or environment containing the
#'   variables in `formula`.
#' @param s an optional \eqn{q \times k}{q x k} 0/1 matrix that permits excluding
#'   some model coefficients (\eqn{q} = number of columns of the model matrix,
#'   \eqn{k} = number of basis functions in `formula.p`).
#' @param nlambda number of values of \eqn{\lambda} (default 100).
#' @param lambda.min.ratio smallest value of \eqn{\lambda}, as a fraction of
#'   \eqn{\lambda_{max}}. The default is `1e-4` if the number of observations is
#'   larger than the number of columns of the model matrix, and `0.01` otherwise.
#' @param lambda an optional user-supplied decreasing sequence of values of
#'   \eqn{\lambda} (on the internal scale, see Details). It overrides `nlambda`
#'   and `lambda.min.ratio`.
#' @param tol convergence tolerance (default `1e-6`).
#' @param maxit maximum number of Newton iterations for each \eqn{\lambda}
#'   (default 100).
#' @param penalty either `"lasso"` (default), the \eqn{\ell_1} penalty
#'   \eqn{\lambda\sum_{j,h} w_j|\theta_{jh}|} of Sottile et al. (2020, eq. 3.4),
#'   or `"group"`, the group lasso penalty
#'   \eqn{\lambda\sum_j w_j\sqrt{k_j}\,\|\theta_{j\cdot}\|_2} that selects
#'   covariates (all the coefficients \eqn{\theta_{j1},\dots,\theta_{jk}} of the
#'   \eqn{j}-th covariate are either zero or nonzero). Groups are the columns of
#'   the model matrix (each dummy variable of a factor is a separate group), and
#'   \eqn{k_j} is the number of free coefficients of column \eqn{j}. Ignored if
#'   `method = "legacy"`.
#' @param penalty.factor optional vector of non-negative penalty weights
#'   \eqn{w_j}, one for each column of the model matrix except the intercept
#'   (which is never penalized). A zero weight leaves the covariate unpenalized
#'   (it is always in the model). The default is 1 for all covariates. Different
#'   weights can be used, for example, to compute an adaptive lasso.
#' @param method `"pnewton"` (default) for the proximal Newton algorithm
#'   described in Details, or `"legacy"` for the algorithm of 'qrcmNP' 0.2.x,
#'   kept for reproducibility of published results.
#' @param verbose logical: if `TRUE`, the progress along the path is printed.
#' @param ... not used (an error is given for the removed argument `display`).
#'
#' @details
#' Under the quantile regression coefficients model (QRCM) of Frumento and
#' Bottai (2016), the conditional quantile function is
#' \deqn{Q(p \mid x, \theta) = x^\top \theta\, b(p),}{Q(p | x, theta) = x' theta b(p),}
#' where \eqn{b(p)} is a set of \eqn{k} known functions of \eqn{p} (defined
#' by `formula.p`) and \eqn{\theta} is a \eqn{q \times k}{q x k} matrix of
#' parameters. The penalized estimator minimizes
#' \deqn{L_{PEN}(\theta) = \int_0^1 \sum_i w_i \rho_p\{y_i - Q(p \mid x_i, \theta)\}\,dp + \lambda\,\mathrm{pen}(\theta),}{%
#' L(theta) + lambda * pen(theta),}
#' where \eqn{\rho_p} is the check loss of quantile regression. As suggested by
#' Sottile et al. (2020), \eqn{x} and \eqn{b(p)} are standardized, and the
#' response is mapped to the interval \eqn{[0, 10]}; the values of `lambda` refer
#' to this internal scale. The estimated coefficients are returned on the
#' original scale.
#'
#' With `method = "pnewton"`, for each \eqn{\lambda} the integrated loss is
#' replaced by its second-order expansion (gradient and Hessian are available in
#' closed form) and the resulting penalized quadratic problem is solved by
#' coordinate descent (`"lasso"`) or block coordinate descent (`"group"`); a
#' backtracking line search guarantees the decrease of the penalized objective
#' (proximal Newton method). The path starts from
#' \eqn{\lambda_{max}}, the smallest value at which all penalized coefficients
#' are zero, and uses warm starts, sequential strong rules and a final check of
#' the Karush-Kuhn-Tucker conditions. With censored or truncated data, the
#' estimating equations of Frumento and Bottai (2017) replace the gradient of the
#' loss; since they are not the gradient of a loss function, the symmetric part
#' of their Jacobian is used in the quadratic approximation, and the iterate with
#' the smallest violation of the KKT conditions is returned (convergence is
#' declared when the average violation is below `5e-4`).
#'
#' The intercept (the first row of \eqn{\theta}) is never penalized. Use
#' [gof.piqr()] to select \eqn{\lambda} by AIC or BIC, [cv_piqr()] for
#' cross-validation, and [refit_piqr()] to re-estimate the selected model without
#' penalty and obtain standard errors (Sottile et al., 2020, Section 3.2).
#'
#' @return An object of class `"piqr"`, a list with items:
#' \item{call}{the matched call.}
#' \item{lambda}{the sequence of values of \eqn{\lambda}.}
#' \item{coefficients}{a list (one element for each \eqn{\lambda}) of
#'   \eqn{q \times k}{q x k} matrices of estimated coefficients.}
#' \item{covar}{a list of covariance matrices (only with `method = "legacy"`,
#'   otherwise `NULL`; use [refit_piqr()] for inference).}
#' \item{obj.function}{the minimized (unpenalized) integrated loss, on the
#'   original scale, for each \eqn{\lambda}. For censored and truncated data
#'   this is the loss function of 'qrcm', which is not the function minimized.}
#' \item{df}{number of covariates (including the intercept) with nonzero
#'   coefficients, for each \eqn{\lambda}.}
#' \item{npar}{number of nonzero coefficients, for each \eqn{\lambda}.}
#' \item{dl}{matrix of the \eqn{\ell_1}-norms of the rows of the (unpenalized)
#'   gradient, one column for each \eqn{\lambda}.}
#' \item{seqS}{a list of 0/1 matrices indicating the nonzero coefficients.}
#' \item{penalty, penalty.factor, method}{as in the call.}
#' \item{converged}{logical vector, convergence at each \eqn{\lambda}.}
#' \item{internal}{a list of objects used by the methods.}
#'
#' @references
#' Sottile, G., Frumento, P., Chiodi, M. and Bottai, M. (2020). A penalized
#' approach to covariate selection through quantile regression coefficient
#' models. *Statistical Modelling*, 20(4), 369-385. \doi{10.1177/1471082X19825523}
#'
#' Frumento, P. and Bottai, M. (2016). Parametric modeling of quantile
#' regression coefficient functions. *Biometrics*, 72(1), 74-84.
#' \doi{10.1111/biom.12410}
#'
#' Frumento, P. and Bottai, M. (2017). Parametric modeling of quantile
#' regression coefficient functions with censored and truncated data.
#' *Biometrics*, 73(4), 1179-1188. \doi{10.1111/biom.12675}
#'
#' Friedman, J., Hastie, T. and Tibshirani, R. (2010). Regularization paths for
#' generalized linear models via coordinate descent. *Journal of Statistical
#' Software*, 33(1), 1-22. \doi{10.18637/jss.v033.i01}
#'
#' Lee, J. D., Sun, Y. and Saunders, M. A. (2014). Proximal Newton-type methods
#' for minimizing composite functions. *SIAM Journal on Optimization*, 24(3),
#' 1420-1443. \doi{10.1137/130921428}
#'
#' @seealso [gof.piqr()], [cv_piqr()], [refit_piqr()], [summary.piqr()],
#'   [plot.piqr()], [predict.piqr()].
#'
#' @examples
#' set.seed(1234)
#' n <- 300
#' x1 <- rexp(n)
#' x2 <- runif(n, 0, 5)
#' x3 <- rnorm(n)                      # irrelevant covariate
#' # true model: beta0(p) = 1 + qnorm(p), beta1(p) = 1 + 2 p, beta2(p) = 1
#' u <- runif(n)
#' y <- 1 + qnorm(u) + (1 + 2 * u) * x1 + x2
#' dat <- data.frame(y, x1, x2, x3)
#'
#' fit <- piqr(y ~ x1 + x2 + x3, formula.p = ~ slp(p, 2), data = dat, nlambda = 30)
#' fit
#' best <- gof.piqr(fit, method = "BIC", plot = FALSE)
#' coef(fit, pos.lambda = best$posMinLambda)
#'
#' # group lasso: covariates enter or leave the model as a whole
#' fitg <- piqr(y ~ x1 + x2 + x3, formula.p = ~ slp(p, 2), data = dat,
#'              nlambda = 30, penalty = "group")
#' plot(fitg, xvar = "lambda")
#'
#' # unpenalized re-estimation of the selected model (with standard errors)
#' \donttest{
#' m <- refit_piqr(fit, pos.lambda = best$posMinLambda)
#' summary(m)
#' }
#' @export
piqr <- function(formula, formula.p = ~ slp(p, 3), weights, data, s, nlambda = 100,
                 lambda.min.ratio = ifelse(nobs < nvars, 0.01, 1e-04), lambda,
                 tol = 1e-6, maxit = 100, penalty = c("lasso", "group"),
                 penalty.factor, method = c("pnewton", "legacy"), verbose = FALSE, ...) {

  cl <- match.call()
  penalty <- match.arg(penalty)
  method <- match.arg(method)
  dots <- names(list(...))
  if ("display" %in% dots)
    stop("the argument 'display' was removed in qrcmNP 1.0.0: please use 'verbose'", call. = FALSE)
  if (length(dots)) stop("unused argument(s): ", paste(dots, collapse = ", "), call. = FALSE)

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "weights", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  A <- check_in(mf, formula.p, if (missing(s)) NULL else s)
  nobs <- nrow(A$V$X)
  nvars <- ncol(A$V$X)
  force(lambda.min.ratio)
  if (!is.numeric(nlambda) || nlambda < 1) stop("'nlambda' must be a positive integer", call. = FALSE)
  if (!missing(lambda)) {
    if (any(lambda < 0)) stop("'lambda' must be non-negative", call. = FALSE)
    lambda <- sort(lambda, decreasing = TRUE)
  } else lambda <- NULL

  pf <- make_penalty_factor(if (missing(penalty.factor)) NULL else penalty.factor, A)

  fit <- if (method == "pnewton") {
    piqr_path(A, penalty = penalty, pf = pf, nlambda = nlambda,
              lambda.min.ratio = lambda.min.ratio, lambda = lambda,
              tol = tol, maxit = maxit, verbose = verbose)
  } else {
    if (!missing(penalty.factor)) warning("'penalty.factor' is ignored with method = \"legacy\"", call. = FALSE)
    piqr_legacy(A, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, lambda = lambda,
                tol = tol, maxit = maxit, verbose = verbose)
  }

  mf <- A$mf
  S <- list(B = A$stats.B, X = A$stats.X, y = A$stats.y)
  attributes(A$bfun) <- c(attributes(A$bfun), S$B)
  attr(mf, "assign") <- S$X$assign
  attr(mf, "stats") <- S
  attr(mf, "all.vars") <- A$V
  attr(mf, "bfun") <- A$bfun
  attr(mf, "type") <- A$type
  attr(mf, "contrasts") <- attr(A$U$X, "contrasts")

  nl <- length(fit$lambda)
  labs <- paste0("lambda", seq_len(nl))
  cn <- list(S$X$coef.names, S$B$coef.names)
  coefficients <- lapply(fit$theta, function(th) { dimnames(th) <- cn; th })
  # selected coefficients: nonzero on the internal (standardised) scale, where the
  # penalty acts (the back-transformation can fill the constant column of b(p))
  seqS <- lapply(fit$theta.int, function(th) { s0 <- 1 * (th != 0); dimnames(s0) <- cn; s0 })
  names(coefficients) <- names(seqS) <- labs
  if (!is.null(fit$covar)) names(fit$covar) <- labs
  dl <- fit$dl; dimnames(dl) <- list(S$X$coef.names, labs)

  out <- list(call = cl, lambda = fit$lambda, coefficients = coefficients, covar = fit$covar,
              obj.function = fit$obj.function,
              df = vapply(seqS, function(x) sum(rowSums(x) > 0), numeric(1)),
              npar = vapply(seqS, sum, numeric(1)),
              dl = dl, seqS = seqS, penalty = if (method == "legacy") "legacy" else penalty,
              penalty.factor = pf, method = method,
              converged = fit$converged)
  names(out$df) <- names(out$npar) <- labs
  out$internal <- list(mf = mf, formula = formula, formula.p = formula.p,
                       nlambda = nl, X = A$U$X, y = A$U$y, k = ncol(A$s),
                       intercept = S$X$intercept, type = A$type, s = A$s, pf = pf)
  attr(out$internal$y, "sd") <- stats::sd(A$U$y)
  class(out) <- "piqr"
  out
}

#' Penalty factor for each column of the model matrix (0 for the intercept).
#' @noRd
make_penalty_factor <- function(penalty.factor, A) {
  q <- ncol(A$V$X)
  intX <- A$stats.X$intercept
  cX <- A$stats.X$const
  nonint <- setdiff(seq_len(q), if (intX) cX else integer(0))
  pf <- rep.int(1, q)
  if (!is.null(penalty.factor)) {
    if (!is.numeric(penalty.factor) || any(penalty.factor < 0) || anyNA(penalty.factor))
      stop("'penalty.factor' must be a vector of non-negative numbers", call. = FALSE)
    if (length(penalty.factor) == 1) penalty.factor <- rep.int(penalty.factor, length(nonint))
    if (length(penalty.factor) != length(nonint))
      stop(sprintf("'penalty.factor' must have length %d (one weight for each column of the model matrix, excluding the intercept)",
                   length(nonint)), call. = FALSE)
    pf[nonint] <- penalty.factor
  }
  if (intX) pf[cX] <- 0
  names(pf) <- A$stats.X$coef.names
  pf
}

#' Regularisation path with the proximal Newton solver.
#' @noRd
piqr_path <- function(A, penalty, pf, nlambda, lambda.min.ratio, lambda, tol, maxit, verbose) {

  V <- A$V; bfun <- A$internal.bfun; type <- A$type
  free <- (A$s == 1)
  q <- nrow(free); k <- ncol(free)
  ee <- get_ee(type)
  S <- list(B = A$stats.B, X = A$stats.X, y = A$stats.y)

  pen <- list(type = penalty, free = free,
              W = matrix(pf, q, k) * free,
              gw = pf * sqrt(rowSums(free)))
  unpen <- free & (matrix(pf, q, k) == 0)
  solver <- pnewton

  # null model: only unpenalized coefficients
  theta <- matrix(0, q, k)
  if (any(unpen)) {
    theta <- start_theta(V, bfun, 1 * unpen, type)
    nf <- solver(theta, V, bfun, type, lambda = 0, pen = pen, active = unpen, tol = tol,
                 maxit = max(maxit, 10 + 10 * sum(unpen)))
    theta <- nf$theta; g <- nf$g
  } else {
    ps <- pstar(theta, V, bfun, type)
    g <- ee(theta, V$y, V$z, V$d, V$X, V$Xw, bfun, ps$y, ps$z, J = FALSE)$g
    nf <- list(theta = theta, ps = ps, g = g, converged = TRUE)
  }
  # metric used when J = 0 (theta = 0 without unpenalized coefficients)
  H0 <- NULL
  if (!any(unpen)) {
    th0 <- start_theta(V, bfun, 1 * free, type)
    ps0 <- pstar(th0, V, bfun, type)
    J0 <- ee(th0, V$y, V$z, V$d, V$X, V$Xw, bfun, ps0$y, ps0$z, J = TRUE)$J
    H0 <- (J0 + t(J0)) / 2
  }

  crit <- function(g) {
    # dual norm of the gradient, for each element (lasso) or row (group)
    G <- matrix(g, q, k) * free
    if (penalty == "lasso") abs(G) / ifelse(pen$W > 0, pen$W, Inf)
    else sqrt(rowSums(G^2)) / ifelse(pen$gw > 0, pen$gw, Inf)
  }
  lmax <- max(c(crit(g), 0))
  if (lmax == 0) lmax <- 1
  user.lambda <- !is.null(lambda)
  if (!user.lambda) {
    lambda <- if (nlambda == 1) lmax
              else exp(seq(log(lmax), log(lambda.min.ratio * lmax), length.out = nlambda))
  }
  nl <- length(lambda)

  thetas <- thetas.int <- vector("list", nl)
  set_active <- function(theta) {
    if (penalty == "group") unpen | (matrix(rowSums(theta != 0) > 0, q, k) & free)
    else unpen | (theta != 0 & free)
  }
  obj <- numeric(nl); conv <- logical(nl)
  dl <- matrix(0, q, nl)
  lprev <- lmax
  active <- set_active(theta)

  for (t in seq_len(nl)) {
    lam <- lambda[t]
    cr <- crit(g)
    strong <- if (penalty == "lasso") (cr >= 2 * lam - lprev) & free
              else matrix(cr >= 2 * lam - lprev, q, k) & free
    act <- active | strong
    if (t == 1 && !user.lambda) act <- unpen   # lambda[1] = lambda_max: the null model
    repeat {
      if (t == 1 && !user.lambda) { fit <- nf; break }
      fit <- solver(theta, V, bfun, type, lambda = lam, pen = pen, active = act,
                    tol = tol, maxit = maxit, H0 = H0)
      theta <- fit$theta; g <- fit$g
      # KKT check on the excluded coefficients
      cr <- crit(g)
      viol <- if (penalty == "lasso") (cr > lam * (1 + 1e-3)) & free & !act
              else matrix(cr > lam * (1 + 1e-3), q, k) & free & !act
      if (!any(viol)) break
      act <- act | viol
    }
    active <- set_active(theta)
    lprev <- lam
    thetas.int[[t]] <- theta
    thetas[[t]] <- check_out(theta, S)$theta
    obj[t] <- loss_internal(theta, V, bfun, type, fit$ps) * (S$y$M - S$y$m) / 10
    conv[t] <- fit$converged
    dl[, t] <- rowSums(abs(matrix(g, q, k) * free))
    if (verbose) {
      message(sprintf("lambda %3d: %10.5f   obj.function = %10.4f   df = %d   npar = %d",
                      t, lam, obj[t], sum(rowSums(theta != 0) > 0), sum(theta != 0)))
    }
  }
  if (!all(conv)) warning(sprintf("the algorithm did not converge for %d value(s) of lambda", sum(!conv)),
                          call. = FALSE)
  list(lambda = lambda, theta = thetas, theta.int = thetas.int, covar = NULL, obj.function = obj,
       dl = dl, converged = conv)
}
