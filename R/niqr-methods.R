# Methods and auxiliary functions for "niqr" objects.

#' Model list used by the engine, from a fitted object (optionally new X, y).
#' @noRd
niqr_M <- function(object, X = object$X, y = object$y) {
  list(X = X, y = y, weights = rep.int(1, nrow(X)), fun = object$fun, fun2 = object$fun2,
       h = object$control$h, grid = niqr_grid(object$control$n.p, object$control$eps), n.bisec = 30)
}

#' Model matrix (and response) for new data.
#' @noRd
niqr_newdata <- function(object, newdata, response = FALSE) {
  if (missing(newdata) || is.null(newdata)) return(list(X = object$X, y = object$y))
  if (!is.data.frame(newdata) && !is.list(newdata))
    stop("'newdata' must be a data frame", call. = FALSE)
  mt <- object$terms
  if (!response) mt <- delete.response(mt)
  mf <- model.frame(mt, data = as.data.frame(newdata), xlev = object$xlevels, na.action = stats::na.exclude)
  list(X = model.matrix(mt, mf, contrasts.arg = object$contrasts),
       y = if (response) model.response(mf, "numeric") else NULL,
       na.action = attr(mf, "na.action"), nrow = nrow(as.data.frame(newdata)))
}

#' Numerical Jacobian of a vector function of theta.
#' @noRd
num_jac <- function(f, theta, h = .Machine$double.eps^(1 / 3)) {
  f0 <- f(theta)
  J <- matrix(NA_real_, length(f0), length(theta))
  for (k in seq_along(theta)) {
    hk <- h * max(1, abs(theta[k]))
    tp <- tm <- theta; tp[k] <- tp[k] + hk; tm[k] <- tm[k] - hk
    J[, k] <- (f(tp) - f(tm)) / (2 * hk)
  }
  J
}

#' @export
print.niqr <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Coefficients:\n")
  print.default(format(rbind(Estimate = x$coefficients, `Std. Error` = x$se), digits = digits),
                print.gap = 2L, quote = FALSE)
  cat("\nMinimized loss function:", format(x$obj.function, digits = digits + 3),
      "   n. of observations:", length(x$y), "\n\n")
  invisible(x)
}

#' @export
coef.niqr <- function(object, ...) object$coefficients

#' @export
vcov.niqr <- function(object, ...) object$covar

#' @export
nobs.niqr <- function(object, ...) length(object$y)

#' @export
terms.niqr <- function(x, ...) x$terms

#' @export
model.matrix.niqr <- function(object, ...) object$X

#' Summarize a Nonlinear Parametric Quantile Model
#'
#' @param object an object of class `"niqr"`.
#' @param p an optional vector of quantile orders: if supplied, the summary
#'   reports the coefficient functions \eqn{\beta(p \mid \hat\theta)}{beta(p | theta)}
#'   (the columns of `fun`) at each `p`, with delta-method standard errors.
#' @param x an object of class `"summary.niqr"`.
#' @param digits number of significant digits.
#' @param ... for future methods.
#' @return An object of class `"summary.niqr"`. Without `p`: the table of
#'   estimated parameters with standard errors, z-values and p-values, the
#'   minimized loss and the information criteria of [gof.niqr()]. With `p`: a
#'   list with one table for each value of `p`.
#' @seealso [niqr()], [gof.niqr()]
#' @export
summary.niqr <- function(object, p, ...) {
  if (missing(p)) {
    est <- object$coefficients; se <- object$se
    z <- est / se
    coe <- cbind(Estimate = est, `Std. Error` = se, `z value` = z, `Pr(>|z|)` = 2 * pnorm(-abs(z)))
    out <- list(call = object$call, coefficients = coe, obj.function = object$obj.function,
                n = length(object$y), n.it = object$n.it, converged = object$converged,
                gof = gof.niqr(object), se.type = object$control$se,
                cluster = !is.null(object$cluster))
  } else {
    if (any(p <= 0 | p >= 1)) stop("0 < p < 1 is required", call. = FALSE)
    b <- predict.niqr(object, type = "beta", p = p, se = TRUE)
    out <- lapply(seq_along(p), function(i) {
      est <- vapply(b, function(bb) bb$beta[i], numeric(1))
      se <- vapply(b, function(bb) if (is.null(bb$se)) NA_real_ else bb$se[i], numeric(1))
      z <- est / se
      cbind(Estimate = est, `Std. Error` = se, `z value` = z, `Pr(>|z|)` = 2 * pnorm(-abs(z)))
    })
    names(out) <- paste("p =", p)
    out <- list(call = object$call, beta = out)
  }
  class(out) <- "summary.niqr"
  out
}

#' @rdname summary.niqr
#' @export
print.summary.niqr <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if (!is.null(x$coefficients)) {
    cat("n. of observations:", x$n, "   n. of iterations:", x$n.it,
        if (!x$converged) "  (not converged)", "\n")
    cat("Standard errors:", switch(x$se.type, sandwich = "sandwich", boot = "bootstrap", none = "none"),
        if (x$cluster) "(clustered)", "\n\n")
    printCoefmat(x$coefficients, digits = digits, signif.stars = TRUE, has.Pvalue = TRUE, P.values = TRUE)
    cat("\nMinimized loss function:", format(x$obj.function, digits = digits + 3), "\n")
    g <- x$gof
    cat("AIC =", format(g$AIC, digits = digits + 2), "  BIC =", format(g$BIC, digits = digits + 2),
        "  AIC* =", format(g$AIC.log, digits = digits + 2), "  BIC* =", format(g$BIC.log, digits = digits + 2),
        "\n\n")
  } else {
    for (j in seq_along(x$beta)) {
      cat(names(x$beta)[j], "\n")
      printCoefmat(x$beta[[j]], digits = digits, signif.stars = TRUE, has.Pvalue = TRUE, P.values = TRUE)
      cat("\n")
    }
  }
  invisible(x)
}

#' Prediction after Nonlinear Parametric Quantile Models
#'
#' @param object an object of class `"niqr"`.
#' @param type the type of prediction: `"beta"` (the coefficient functions,
#'   i.e. the columns of `fun(theta, p)`), `"CDF"` and `"PDF"` (fitted
#'   conditional distribution and density at the observed response), `"QF"`
#'   (conditional quantile function) or `"sim"` (simulated responses).
#' @param newdata an optional data frame (for `type = "CDF"` and `"PDF"`
#'   it must contain the response).
#' @param p a vector of quantile orders (for `type = "beta"` and `"QF"`).
#' @param se logical: compute delta-method standard errors (`"beta"` and `"QF"`)?
#' @param ... not used.
#' @return
#' * `type = "beta"`: a list with one data frame for each column of `fun`, with
#'   columns `p`, `beta`, and (if `se = TRUE`) `se`, `low`, `up` (95% pointwise
#'   confidence limits).
#' * `type = "CDF"` or `"PDF"`: a vector.
#' * `type = "QF"`: a data frame with one column for each `p` or, if `se = TRUE`,
#'   a list with items `fit` and `se.fit`.
#' * `type = "sim"`: a vector of simulated responses.
#' @seealso [niqr()]
#' @export
predict.niqr <- function(object, type = c("beta", "CDF", "PDF", "QF", "sim"), newdata, p,
                         se = TRUE, ...) {
  type <- match.arg(type)
  theta <- object$coefficients
  V <- object$covar
  if (se && all(is.na(V))) se <- FALSE

  if (type == "beta") {
    if (missing(p)) p <- seq.int(0.01, 0.99, 0.01)
    if (any(p <= 0 | p >= 1)) stop("0 < p < 1 is required", call. = FALSE)
    B <- object$fun(theta, p)
    if (!is.matrix(B)) B <- matrix(B, nrow = length(p))
    nm <- colnames(B)
    if (is.null(nm)) nm <- if (is.null(object$fun2) && ncol(B) == ncol(object$X)) colnames(object$X)
                           else paste0("beta", seq_len(ncol(B)) - 1)
    out <- lapply(seq_len(ncol(B)), function(j) {
      d <- data.frame(p = p, beta = B[, j])
      if (se) {
        J <- num_jac(function(th) { b <- object$fun(th, p); matrix(b, nrow = length(p))[, j] }, theta,
                     object$control$h)
        d$se <- sqrt(pmax(rowSums((J %*% V) * J), 0))
        d$low <- d$beta - 1.96 * d$se
        d$up <- d$beta + 1.96 * d$se
      }
      d
    })
    names(out) <- nm
    return(out)
  }

  nd <- niqr_newdata(object, if (missing(newdata)) NULL else newdata, response = type %in% c("CDF", "PDF"))
  M <- niqr_M(object, nd$X, if (is.null(nd$y)) rep(0, nrow(nd$X)) else nd$y)

  pad <- function(v) {
    na <- nd$na.action
    if (is.null(na)) return(v)
    out <- if (is.matrix(v) || is.data.frame(v)) {
      o <- matrix(NA_real_, nd$nrow, ncol(v)); o[-na, ] <- as.matrix(v); o
    } else { o <- rep(NA_real_, nd$nrow); o[-na] <- v; o }
    out
  }
  if (type %in% c("CDF", "PDF")) {
    if (missing(newdata)) return(if (type == "CDF") object$CDF else object$PDF)
    ev <- niqr_eval(theta, M, hessian = TRUE)
    return(pad(if (type == "CDF") c(ev$ps) else ev$pdf))
  }
  if (type == "QF") {
    if (missing(p)) stop("please indicate the value(s) of 'p'", call. = FALSE)
    if (any(p <= 0 | p >= 1)) stop("0 < p < 1 is required", call. = FALSE)
    fit <- sef <- matrix(NA_real_, nrow(nd$X), length(p))
    for (j in seq_along(p)) {
      qd <- niqr_qd(theta, p[j], M, deriv = se)
      fit[, j] <- qd$Q
      if (se) sef[, j] <- sqrt(pmax(rowSums((qd$D %*% V) * qd$D), 0))
    }
    fit <- as.data.frame(pad(fit)); colnames(fit) <- paste0("p", p)
    if (se) { sef <- as.data.frame(pad(sef)); colnames(sef) <- colnames(fit); return(list(fit = fit, se.fit = sef)) }
    return(fit)
  }
  # sim
  pad(niqr_qd(theta, runif(nrow(nd$X)), M)$Q)
}

#' Plot a Nonlinear Parametric Quantile Model
#'
#' Plots the estimated coefficient functions (the columns of `fun(theta, p)`)
#' with pointwise confidence bands, and a uniform QQ-plot of the fitted CDF.
#'
#' @param x an object of class `"niqr"`.
#' @param conf.int logical: draw 95% confidence bands?
#' @param which an optional vector of indices of the coefficient functions to
#'   plot (use `length(beta) + 1` for the QQ-plot). If missing and
#'   `ask = TRUE`, a menu is shown.
#' @param ask logical: choose interactively what to plot? Default `interactive()`.
#' @param polygon logical: draw the confidence bands as shaded areas?
#' @param ... additional graphical parameters.
#' @return No return value, called for its side effects.
#' @seealso [niqr()]
#' @export
plot.niqr <- function(x, conf.int = TRUE, which = NULL, ask = interactive(), polygon = TRUE, ...) {
  L <- list(...)
  if (is.null(L$xlim)) L$xlim <- c(0.01, 0.99)
  if (is.null(L$lwd)) L$lwd <- 2
  if (is.null(L$col)) L$col <- "black"
  if (is.null(L$xlab)) L$xlab <- "p"
  p <- seq.int(max(0.001, L$xlim[1]), min(0.999, L$xlim[2]), length.out = 100)
  u <- predict.niqr(x, type = "beta", p = p, se = conf.int)
  conf.int <- conf.int && !is.null(u[[1]]$se)
  labels <- if (is.null(L$labels)) names(u) else L$labels
  q <- length(u)
  if (is.null(L$ylab)) L$ylab <- rep("beta(p)", q)
  one <- function(j) {
    ylim <- L$ylim
    if (is.null(ylim)) ylim <- if (conf.int) range(u[[j]]$low, u[[j]]$up) else range(u[[j]]$beta)
    plot(p, u[[j]]$beta, type = "l", lwd = L$lwd, col = L$col, xlab = L$xlab, ylab = L$ylab[j],
         main = labels[j], xlim = L$xlim, ylim = ylim)
    if (conf.int) {
      if (polygon) polygon(c(p, rev(p)), c(u[[j]]$low, rev(u[[j]]$up)),
                           col = adjustcolor(L$col, alpha.f = 0.25), border = NA)
      else { lines(p, u[[j]]$low, lty = 2, col = L$col); lines(p, u[[j]]$up, lty = 2, col = L$col) }
    }
  }
  qq <- function() {
    cdf <- sort(x$CDF)
    plot(ppoints(length(cdf)), cdf, pch = 20, cex = 0.5, xlim = c(0, 1), ylim = c(0, 1),
         xlab = "U(0,1) quantiles", ylab = "fitted CDF quantiles")
    abline(0, 1)
  }
  if (!is.null(which) || !ask) {
    if (is.null(which)) which <- seq_len(q)
    for (j in which) if (j <= q) one(j) else qq()
  } else {
    pick <- 1
    while (pick > 0 && pick <= q + 1) {
      pick <- menu(c(labels, "qqplot"), title = "Make a plot selection (or 0 to exit):\n")
      if (pick > 0 && pick <= q) one(pick) else if (pick == q + 1) qq()
    }
  }
  invisible(u)
}

#' Information Criteria for Nonlinear Parametric Quantile Models
#'
#' Computes the goodness-of-fit criteria of Bottai and Cilluffo (2020,
#' Section 3) for one or more [niqr()] fits.
#'
#' @param object an object of class `"niqr"`.
#' @param ... further objects of class `"niqr"` fitted to the same data.
#' @details With \eqn{L_n} the minimized (average) loss, \eqn{k} the number of
#'   parameters and \eqn{n} the sample size, the criteria are
#'   \eqn{AIC^* = \log L_n + k/n} and \eqn{BIC^* = \log L_n + k\log(n)/(2n)}
#'   (`AIC.log`, `BIC.log`), and \eqn{AIC = nL_n + 2k}, \eqn{BIC = nL_n + k\log n}.
#'   Smaller values indicate a better fit. Bottai and Cilluffo (2020) write the
#'   last two criteria in terms of \eqn{L_n}; here the total loss \eqn{nL_n} is
#'   used, otherwise the penalty would dominate. Note that, unlike the log
#'   versions, `AIC` and `BIC` depend on the scale of the response; the log
#'   versions are recommended.
#' @return A data frame with columns `npar`, `loss`, `AIC`, `BIC`, `AIC.log`
#'   and `BIC.log` (one row per model).
#' @references Bottai, M. and Cilluffo, G. (2020). Nonlinear parametric
#'   quantile models. *Statistical Methods in Medical Research*, 29(12),
#'   3757-3769. \doi{10.1177/0962280220941159}
#' @seealso [niqr()], [test.fit.niqr()]
#' @export
gof.niqr <- function(object, ...) {
  objs <- c(list(object), list(...))
  if (!all(vapply(objs, inherits, TRUE, what = "niqr"))) stop("all objects must be of class 'niqr'", call. = FALSE)
  out <- do.call(rbind, lapply(objs, function(o) {
    L <- o$obj.function; k <- length(o$coefficients); n <- length(o$y)
    data.frame(npar = k, loss = L, AIC = n * L + 2 * k, BIC = n * L + k * log(n),
               AIC.log = log(L) + k / n, BIC.log = log(L) + k * log(n) / (2 * n))
  }))
  if (length(objs) > 1) {
    nm <- vapply(as.list(match.call())[-1], function(e) paste(deparse(e), collapse = ""), "")
    rownames(out) <- nm
    if (length(unique(vapply(objs, function(o) length(o$y), 1))) > 1)
      warning("the models are not fitted to the same number of observations", call. = FALSE)
  }
  out
}

#' Goodness-of-Fit Test for Nonlinear Parametric Quantile Models
#'
#' Under a correctly specified model, the fitted values of the conditional
#' CDF, \eqn{F(y_i \mid x_i, \hat\theta)}{F(y_i | x_i, theta)}, are approximately
#' uniformly distributed. This function computes the Kolmogorov-Smirnov and
#' Cramer-von Mises statistics of their empirical distribution, and obtains
#' p-values by parametric bootstrap, as in [qrcm::test.fit()].
#'
#' @param object an object of class `"niqr"`.
#' @param R number of bootstrap replications (default 100).
#' @param trace logical: print the progress?
#' @param ... not used.
#' @return A matrix with the test statistics and their bootstrap p-values.
#' @seealso [niqr()], [gof.niqr()]
#' @export
test.fit.niqr <- function(object, R = 100, trace = FALSE, ...) {
  unif_stat <- function(Fy) {
    n <- length(Fy)
    Fy <- sort(Fy)
    hat <- seq_len(n) / n
    ks <- max(hat - Fy, Fy - (seq_len(n) - 1) / n)
    Fy2 <- c(0, Fy, 1); hat2 <- c(0, hat, 1)
    U <- (hat2 - Fy2)^2
    h <- diff(Fy2)
    cvm <- sum((U[-1] + U[-length(U)]) * h / 2)
    c(ks = ks, cvm = cvm)
  }
  n <- length(object$y)
  t0 <- unif_stat(object$CDF)
  ctrl <- object$control; ctrl$se <- "none"; ctrl$verbose <- FALSE
  Tb <- matrix(NA_real_, R, 2)
  M <- niqr_M(object)
  for (b in seq_len(R)) {
    if (trace) message("bootstrap replication ", b, " of ", R)
    id <- sample.int(n, n, replace = TRUE)
    Xb <- object$X[id, , drop = FALSE]
    Mb <- M; Mb$X <- Xb
    yb <- niqr_qd(object$coefficients, runif(n), Mb)$Q
    f <- tryCatch(suppressWarnings(niqr_fit(object$coefficients, yb, Xb, object$weights[id], object$fun,
                                            object$fun2, ctrl, boot = TRUE)), error = function(e) NULL)
    if (!is.null(f) && f$converged) Tb[b, ] <- unif_stat(f$CDF)
  }
  ok <- stats::complete.cases(Tb)
  if (sum(!ok) > 0) warning(sprintf("%d bootstrap replications failed", sum(!ok)), call. = FALSE)
  out <- cbind(statistic = t0 * c(1, n),
               `p-value` = c(mean(Tb[ok, 1] >= t0[1]), mean(Tb[ok, 2] >= t0[2])))
  rownames(out) <- c("Kolmogorov-Smirnov", "Cramer-Von Mises")
  out
}

#' Defunct Functions in qrcmNP
#'
#' These functions and interfaces were removed in 'qrcmNP' 1.0.0:
#' * `testfit.niqr(obj, R)`: use [test.fit()] (`test.fit(object, R)`);
#' * `niqr(fun, fun2, x0, X, y, control)`: use
#'   `niqr(formula, data, fun, fun2, x0)` (see [niqr()]);
#' * the argument `display` of [piqr()]: use `verbose`;
#' * the argument `df.new` of [gof.piqr()]: use `df`;
#' * the control names of 'qrcmNP' 0.2.x in [niqr()]: use [niqr.control()].
#'
#' @param ... ignored.
#' @return No return value: an error is raised.
#' @name qrcmNP-defunct
#' @keywords internal
NULL

#' @rdname qrcmNP-defunct
#' @export
testfit.niqr <- function(...) {
  .Defunct("test.fit", package = "qrcmNP",
           msg = "'testfit.niqr()' was removed in qrcmNP 1.0.0: please use 'test.fit(object, R)'")
}
