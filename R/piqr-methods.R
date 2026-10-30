# Methods and auxiliary functions for "piqr" objects.

#' Index of the selected value of lambda.
#' @noRd
lambda_index <- function(object, pos.lambda, lambda) {
  nl <- length(object$lambda)
  if (!missing(pos.lambda) && !is.null(pos.lambda)) {
    l <- as.integer(pos.lambda)[1]
    if (is.na(l) || l < 1 || l > nl) stop(sprintf("'pos.lambda' must be an integer between 1 and %d", nl), call. = FALSE)
    return(l)
  }
  if (!missing(lambda) && !is.null(lambda)) return(which.min(abs(object$lambda - lambda[1])))
  if (nl == 1) return(1L)
  stop("please indicate the position of a lambda value with 'pos.lambda' ",
       "(see 'gof.piqr()' to select it by AIC or BIC)", call. = FALSE)
}

#' Model list (coefficients, covar, mf) at a given lambda, as used by the engine.
#' @noRd
piqr_model <- function(object, l) {
  list(coefficients = object$coefficients[[l]],
       covar = if (is.null(object$covar)) NULL else object$covar[[l]],
       mf = object$internal$mf)
}

#' @export
print.piqr <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Penalty: ", x$penalty, "   n. of lambda values: ", length(x$lambda), "\n\n", sep = "")
  print.default(cbind(Df = x$df, Npar = x$npar, ObjFunction = signif(x$obj.function, digits),
                      Lambda = signif(x$lambda, digits)), print.gap = 2L, quote = FALSE)
  cat("\n")
  invisible(x)
}

#' Summarize a Penalized QRCM
#'
#' Summary of a [piqr()] fit at a selected value of \eqn{\lambda}.
#'
#' @param object an object of class `"piqr"`.
#' @param pos.lambda the position of the selected value of \eqn{\lambda} in
#'   `object$lambda` (e.g., `gof.piqr(object)$posMinLambda`).
#' @param SE logical: if `TRUE`, standard errors of \eqn{\theta} are reported
#'   (available only with `method = "legacy"`; they ignore the selection step,
#'   see [refit_piqr()]).
#' @param p an optional vector of quantile orders. If supplied, the summary
#'   reports \eqn{\beta(p)} at each \eqn{p}.
#' @param cov logical: if `TRUE` and `p` is supplied, the covariance matrix of
#'   \eqn{\beta(p)} is reported (only with `method = "legacy"`).
#' @param lambda alternatively to `pos.lambda`, a value of \eqn{\lambda}; the
#'   closest value in `object$lambda` is used.
#' @param x an object of class `"summary.piqr"`.
#' @param digits number of significant digits.
#' @param ... for future methods.
#'
#' @return An object of class `"summary.piqr"`. Without `p`: a list with the
#'   selected `lambda`, the matrix of `coefficients`, their standard errors
#'   `se` (if requested and available), the minimized loss `obj.function`, the
#'   number of observations `n` and of nonzero parameters `free.par`. With
#'   `p`: a list with one element for each value of `p`, containing the
#'   estimated \eqn{\beta(p)} (and standard errors, if available).
#'
#' @seealso [piqr()], [gof.piqr()], [refit_piqr()].
#' @examples
#' set.seed(1)
#' n <- 200
#' x1 <- runif(n); x2 <- runif(n)
#' y <- 1 + qnorm(runif(n)) + 2 * x1
#' fit <- piqr(y ~ x1 + x2, formula.p = ~ slp(p, 2), nlambda = 20)
#' best <- gof.piqr(fit, plot = FALSE)
#' summary(fit, pos.lambda = best$posMinLambda)
#' summary(fit, pos.lambda = best$posMinLambda, p = c(0.1, 0.5, 0.9))
#' @export
summary.piqr <- function(object, pos.lambda, SE = FALSE, p, cov = FALSE, lambda, ...) {
  l <- lambda_index(object, if (missing(pos.lambda)) NULL else pos.lambda,
                    if (missing(lambda)) NULL else lambda)
  model <- piqr_model(object, l)
  if (missing(p)) {
    se <- NULL
    if (SE) {
      if (is.null(model$covar)) {
        message("standard errors are available only with method = \"legacy\"; ",
                "use refit_piqr() for post-selection inference")
      } else {
        se <- matrix(sqrt(pmax(diag(model$covar), 0)), ncol = object$internal$k,
                     dimnames = dimnames(model$coefficients))
      }
    }
    out <- list(lambda = object$lambda[l], pos.lambda = l, coefficients = model$coefficients, se = se,
                covar = model$covar, obj.function = object$obj.function[l],
                n = nrow(object$internal$mf), free.par = object$npar[l], df = object$df[l],
                penalty = object$penalty, call = object$call, SE = !is.null(se))
  } else {
    if (any(p <= 0 | p >= 1)) stop("0 < p < 1 is required", call. = FALSE)
    out <- lapply(p, function(pp) extract_p(model, pp, cov = cov && !is.null(model$covar)))
    names(out) <- paste("p =", p)
    attr(out, "nacoef") <- which(apply(model$coefficients, 1, function(v) all(v == 0)))
    out$call <- object$call
  }
  class(out) <- "summary.piqr"
  out
}

#' @rdname summary.piqr
#' @export
print.summary.piqr <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if (!is.null(x$coefficients)) {
    coe <- x$coefficients
    coe[coe == 0] <- NA
    cat("Penalty:", x$penalty, "\n")
    cat("Selected value of lambda:", format(x$lambda, digits = digits),
        paste0("(position ", x$pos.lambda, ")"), "\n")
    cat("n. of observations:", x$n, "\n")
    cat("n. of selected covariates (incl. intercept):", x$df, "\n")
    cat("n. of nonzero parameters:", x$free.par, "\n\n")
    cat("Coefficients (zero coefficients are shown as NA):\n")
    print.default(format(coe, digits = digits), print.gap = 2L, quote = FALSE)
    cat("\n")
    if (isTRUE(x$SE)) {
      se <- x$se; se[is.na(coe)] <- NA
      cat("Standard errors:\n")
      print.default(format(se, digits = digits), print.gap = 2L, quote = FALSE)
      cat("\n")
    }
    cat("Minimized loss function:", format(x$obj.function, digits = digits + 3), "\n\n")
  } else {
    nacoef <- attr(x, "nacoef")
    for (j in seq_len(length(x) - 1)) {
      cat(names(x)[j], "\n\n")
      coe <- x[[j]]$coef
      coe[nacoef, ] <- NA
      if (ncol(coe) == 4) {
        printCoefmat(coe, digits = digits, signif.stars = TRUE, cs.ind = 1:2, tst.ind = 3,
                     P.values = TRUE, has.Pvalue = TRUE, na.print = "-")
      } else {
        print.default(format(coe, digits = digits), print.gap = 2L, quote = FALSE)
      }
      if (!is.null(x[[j]]$cov)) {
        cat("\nCovariance:\n")
        print.default(format(x[[j]]$cov, digits = digits), print.gap = 2L, quote = FALSE)
      }
      cat("\n")
    }
  }
  invisible(x)
}

#' Extract Coefficients of a Penalized QRCM
#'
#' @param object an object of class `"piqr"`.
#' @param pos.lambda the position of the selected value of \eqn{\lambda}.
#' @param lambda alternatively, a value of \eqn{\lambda} (the closest value in
#'   `object$lambda` is used).
#' @param ... not used.
#' @return The \eqn{q \times k}{q x k} matrix \eqn{\hat\theta} at the selected
#'   \eqn{\lambda}.
#' @seealso [piqr()]
#' @export
coef.piqr <- function(object, pos.lambda, lambda, ...) {
  l <- lambda_index(object, if (missing(pos.lambda)) NULL else pos.lambda,
                    if (missing(lambda)) NULL else lambda)
  object$coefficients[[l]]
}

#' @export
terms.piqr <- function(x, ...) attr(x$internal$mf, "terms")

#' @export
model.matrix.piqr <- function(object, ...) {
  mf <- object$internal$mf
  model.matrix(attr(mf, "terms"), mf)
}

#' @export
nobs.piqr <- function(object, ...) nrow(object$internal$mf)

#' Prediction after Penalized QRCM
#'
#' Predictions from a [piqr()] fit at a selected value of \eqn{\lambda}.
#'
#' @param object an object of class `"piqr"`.
#' @param pos.lambda the position of the selected value of \eqn{\lambda}.
#' @param type the type of prediction: `"beta"` (quantile regression
#'   coefficients \eqn{\beta(p)}), `"CDF"` (fitted conditional distribution and
#'   density functions at the observed response), `"QF"` (conditional quantile
#'   function) or `"sim"` (simulated responses).
#' @param newdata an optional data frame in which to look for the covariates
#'   (and, for `type = "CDF"`, the response). If omitted, the data used to fit
#'   the model are used.
#' @param p a vector of quantile orders (for `type = "beta"` and `type = "QF"`).
#' @param se logical: should standard errors be computed? They are available only
#'   with `method = "legacy"`; see [refit_piqr()] otherwise.
#' @param lambda alternatively to `pos.lambda`, a value of \eqn{\lambda}.
#' @param ... not used.
#' @return
#' * `type = "beta"`: a list with one data frame for each covariate, with
#'   columns `p`, `beta` (and `se`, `low`, `up` if `se = TRUE`).
#' * `type = "CDF"`: a data frame with columns `CDF` and `PDF`.
#' * `type = "QF"`: a data frame of fitted quantiles (one column for each `p`)
#'   or, if `se = TRUE`, a list with items `fit` and `se.fit`.
#' * `type = "sim"`: a vector of simulated responses.
#' @seealso [piqr()]
#' @examples
#' set.seed(1)
#' n <- 200
#' x <- runif(n)
#' y <- qexp(runif(n)) * (1 + x)
#' fit <- piqr(y ~ x, formula.p = ~ slp(p, 3), nlambda = 10)
#' predict(fit, pos.lambda = 10, type = "QF", p = c(0.25, 0.5),
#'         newdata = data.frame(x = c(0, 0.5, 1)))
#' head(predict(fit, pos.lambda = 10, type = "CDF"))
#' @export
predict.piqr <- function(object, pos.lambda, type = c("beta", "CDF", "QF", "sim"), newdata, p,
                         se = !is.null(object$covar), lambda, ...) {
  l <- lambda_index(object, if (missing(pos.lambda)) NULL else pos.lambda,
                    if (missing(lambda)) NULL else lambda)
  type <- match.arg(type)
  model <- piqr_model(object, l)
  if (se && is.null(model$covar)) {
    message("standard errors are available only with method = \"legacy\"; use refit_piqr()")
    se <- FALSE
  }
  theta <- model$coefficients
  bfun <- attr(object$internal$mf, "bfun")

  if (type == "beta") {
    if (missing(p)) p <- seq.int(0.01, 0.99, 0.01)
    if (any(p <= 0 | p >= 1)) stop("0 < p < 1 is required", call. = FALSE)
    return(pred_beta(model, p, se = se))
  }

  mf <- object$internal$mf
  mt <- terms(mf)
  miss <- attr(mf, "na.action")
  fittype <- attr(attr(mf, "type"), "fittype")
  nomiss <- if (is.null(miss)) seq_len(nrow(mf)) else seq_len(nrow(mf) + length(miss))[-miss]
  xlev <- .getXlevels(mt, mf)
  contr <- attr(mf, "contrasts")
  if (!missing(newdata)) {
    newdata <- as.data.frame(newdata)
    if (type == "CDF") {
      yv <- all.vars(mt[[2]])
      yn <- if (fittype == "ctiqr") yv[2] else yv[1]
      if (is.na(match(yn, colnames(newdata))))
        stop("for 'type = \"CDF\"', 'newdata' must contain the response variable", call. = FALSE)
      if (fittype == "ciqr") newdata[, yv[2]] <- 1
      if (fittype == "ctiqr") { newdata[, yv[3]] <- 1; newdata[, yv[1]] <- -Inf }
    } else mt <- delete.response(mt)
    if (any(is.na(match(all.vars(mt), colnames(newdata)))))
      stop("'newdata' must contain all the covariates", call. = FALSE)
    mf <- model.frame(mt, data = newdata, xlev = xlev, na.action = stats::na.exclude)
    miss <- attr(mf, "na.action")
    nomiss <- if (is.null(miss)) seq_len(nrow(mf)) else seq_len(nrow(newdata))[-miss]
  }
  x <- model.matrix(mt, mf, contrasts.arg = contr)
  N <- length(nomiss) + length(miss)

  if (type == "CDF") {
    y <- cbind(model.response(mf))[, 1 + (fittype == "ctiqr")]
    Fy <- p_bisec(theta, y, x, bfun)
    b1 <- apply_bfun(bfun, Fy, "b1fun")
    fy <- 1 / c(rowSums((x %*% theta) * b1))
    if (any(fy < 0)) warning("some PDF values are negative (quantile crossing)", call. = FALSE)
    CDF <- PDF <- rep(NA_real_, N)
    CDF[nomiss] <- Fy; PDF[nomiss] <- fy
    out <- data.frame(CDF = CDF, PDF = PDF)
    rownames(out)[nomiss] <- rownames(mf)
    if (!is.null(miss)) rownames(out)[miss] <- names(miss)
    return(out)
  }
  if (type == "QF") {
    if (missing(p)) stop("please indicate the value(s) of 'p' to compute x'beta(p)", call. = FALSE)
    if (any(p <= 0 | p >= 1)) stop("0 < p < 1 is required", call. = FALSE)
    fit <- se.fit <- matrix(NA_real_, N, length(p))
    for (j in seq_along(p)) {
      fb <- extract_p(model, p[j], cov = se)
      fit[nomiss, j] <- x %*% cbind(fb$coef[, 1])
      if (se) se.fit[nomiss, j] <- sqrt(pmax(rowSums((x %*% fb$cov) * x), 0))
    }
    fit <- as.data.frame(fit); colnames(fit) <- paste0("p", p)
    rownames(fit)[nomiss] <- rownames(mf)
    if (!is.null(miss)) rownames(fit)[miss] <- names(miss)
    if (se) {
      se.fit <- as.data.frame(se.fit); dimnames(se.fit) <- dimnames(fit)
      return(list(fit = fit, se.fit = se.fit))
    }
    return(fit)
  }
  # type == "sim"
  pp <- runif(nrow(x))
  beta <- apply_bfun(bfun, pp, "bfun") %*% t(theta)
  y <- rep(NA_real_, N)
  y[nomiss] <- rowSums(beta * x)
  y
}

#' Plot a Penalized QRCM
#'
#' Coefficient paths, gradient paths and estimated coefficient functions
#' \eqn{\hat\beta(p)} of a [piqr()] fit.
#'
#' @param x an object of class `"piqr"`.
#' @param xvar what to plot: `"lambda"` (coefficient paths against
#'   \eqn{\log\lambda}), `"objective"` (coefficient paths against the minimized
#'   loss), `"grad"` (\eqn{\ell_1}-norm of the gradient of each covariate
#'   against \eqn{\log\lambda}), or `"beta"` (the coefficient functions
#'   \eqn{\hat\beta_j(p)} at `pos.lambda`).
#' @param pos.lambda the position of the selected \eqn{\lambda} (for `xvar = "beta"`).
#' @param label logical: label the curves (paths only)?
#' @param which for `xvar = "beta"`, an optional vector of indices of the
#'   coefficients to plot. If missing and `ask = TRUE`, a menu is shown.
#' @param ask logical: for `xvar = "beta"`, choose interactively which
#'   coefficient to plot? Defaults to `interactive()`.
#' @param polygon logical: for `xvar = "beta"` with standard errors, draw
#'   confidence bands as shaded areas?
#' @param ... additional graphical parameters (`xlab`, `ylab`, `col`, `lwd`,
#'   `xlim`, `ylim`, `main`, ...).
#' @return No return value, called for its side effects.
#' @seealso [piqr()]
#' @examples
#' set.seed(1)
#' n <- 200
#' x1 <- runif(n); x2 <- runif(n); x3 <- runif(n)
#' y <- 1 + qnorm(runif(n)) + 2 * x1 - x2
#' fit <- piqr(y ~ x1 + x2 + x3, formula.p = ~ slp(p, 2), nlambda = 20)
#' oldpar <- par(mfrow = c(1, 2))
#' plot(fit, xvar = "lambda")
#' plot(fit, xvar = "grad")
#' par(oldpar)
#' plot(fit, xvar = "beta", pos.lambda = 20, which = 2)
#' @export
plot.piqr <- function(x, xvar = c("lambda", "objective", "grad", "beta"), pos.lambda,
                      label = FALSE, which = NULL, ask = interactive(), polygon = TRUE, ...) {
  xvar <- match.arg(xvar)
  if (length(x$lambda) == 1 && xvar != "beta") {
    warning("only one value of lambda: xvar = \"beta\" is used", call. = FALSE)
    xvar <- "beta"
  }
  if (xvar == "beta") {
    l <- lambda_index(x, if (missing(pos.lambda)) NULL else pos.lambda)
    return(invisible(plot_beta_piqr(x, l, which = which, ask = ask, polygon = polygon, ...)))
  }
  plot_path_piqr(x, xvar = xvar, label = label, ...)
  invisible(NULL)
}

#' @noRd
plot_path_piqr <- function(x, xvar, label = FALSE, ...) {
  L <- list(...)
  int <- if (x$internal$intercept) attr(x$internal$mf, "stats")$X$const else integer(0)
  if (xvar == "grad") {
    Y <- x$dl[setdiff(seq_len(nrow(x$dl)), int), , drop = FALSE]
    ylab <- "Gradient"
  } else {
    Y <- vapply(x$coefficients, function(th) c(th[setdiff(seq_len(nrow(th)), int), , drop = FALSE]),
                numeric((nrow(x$coefficients[[1]]) - length(int)) * ncol(x$coefficients[[1]])))
    Y <- matrix(Y, ncol = length(x$lambda))
    th1 <- x$coefficients[[1]]
    rn <- outer(rownames(th1)[setdiff(seq_len(nrow(th1)), int)], colnames(th1), paste, sep = ":")
    rownames(Y) <- c(rn)
    Y <- Y[rowSums(abs(Y)) > 0, , drop = FALSE]
    ylab <- "Coefficients"
  }
  if (nrow(Y) == 0) { warning("no plot produced: all coefficients are zero", call. = FALSE); return(invisible()) }
  index <- if (xvar == "objective") x$obj.function else log(x$lambda)
  xlab <- if (xvar == "objective") "Objective function" else expression(log(lambda))
  args <- list(x = index, y = t(Y), type = "l", lty = 1, xlab = xlab, ylab = ylab)
  args[names(L)] <- L
  do.call(matplot, args)
  atdf <- pretty(index)
  prettydf <- approx(x = index, y = x$df, xout = atdf, rule = 2, method = "constant",
                     f = if (xvar == "objective") 1 else 0)$y
  axis(3, at = atdf, labels = prettydf, tcl = NA)
  if (label) {
    xpos <- if (xvar == "objective") max(index) else min(index)
    text(rep(xpos, nrow(Y)), Y[, ncol(Y)], rownames(Y), cex = 0.5, pos = if (xvar == "objective") 4 else 2)
  }
}

#' @noRd
plot_beta_piqr <- function(x, l, which = NULL, ask = interactive(), polygon = TRUE, conf.int = TRUE, ...) {
  L <- list(...)
  if (is.null(L$xlim)) L$xlim <- c(0.01, 0.99)
  if (is.null(L$lwd)) L$lwd <- 2
  if (is.null(L$col)) L$col <- "black"
  if (is.null(L$xlab)) L$xlab <- "p"
  labels <- if (is.null(L$labels)) rownames(x$coefficients[[l]]) else L$labels
  q <- length(labels)
  if (is.null(L$ylab)) L$ylab <- rep("beta(p)", q)
  p <- seq.int(max(0.001, L$xlim[1]), min(0.999, L$xlim[2]), length.out = 100)
  se <- conf.int && !is.null(x$covar)
  u <- predict.piqr(x, pos.lambda = l, p = p, type = "beta", se = se)

  one <- function(j) {
    ylim <- L$ylim
    if (is.null(ylim)) ylim <- if (se) range(u[[j]]$low, u[[j]]$up) else range(u[[j]]$beta)
    plot(p, u[[j]]$beta, type = "l", lwd = L$lwd, col = L$col, xlab = L$xlab, ylab = L$ylab[j],
         main = labels[j], xlim = L$xlim, ylim = ylim)
    if (se) {
      if (polygon) {
        polygon(c(p, rev(p)), c(u[[j]]$low, rev(u[[j]]$up)),
                col = adjustcolor(L$col, alpha.f = 0.25), border = NA)
      } else {
        lines(p, u[[j]]$low, lty = 2, lwd = L$lwd, col = L$col)
        lines(p, u[[j]]$up, lty = 2, lwd = L$lwd, col = L$col)
      }
    }
  }
  qq <- function() {
    cdf <- sort(predict.piqr(x, pos.lambda = l, type = "CDF")$CDF)
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

#' Select the Tuning Parameter of a Penalized QRCM by AIC or BIC
#'
#' Computes an information criterion along the path of a [piqr()] fit and
#' returns the value of \eqn{\lambda} at which it is minimized.
#'
#' @param object an object of class `"piqr"`.
#' @param method the criterion, `"BIC"` (default) or `"AIC"`.
#' @param Cn a positive constant (or a character string with an expression in
#'   `n`, the sample size, and `q`, the number of columns of the model matrix,
#'   e.g. `"log(log(n))"`) multiplying the BIC penalty, see Details.
#' @param plot logical: plot the criterion against \eqn{\log\lambda}?
#' @param logi logical: if `TRUE` (default) the criteria are computed on the
#'   logarithm of the loss (Sottile et al., 2020), otherwise on the loss.
#' @param df how to count the degrees of freedom: `"parameters"`, the number of
#'   nonzero coefficients \eqn{\theta_{jh}} (Sottile et al., 2020, Section 3.1), or
#'   `"covariates"`, the number of covariates with at least one nonzero
#'   coefficient (including the intercept), as in 'qrcmNP' 0.2.x. The default
#'   is `"parameters"`, or `"covariates"` for fits with `method = "legacy"`.
#' @param ... further graphical parameters.
#'
#' @details
#' With `logi = TRUE`, the criteria of Sottile et al. (2020, eq. 3.6-3.7) are
#' \deqn{AIC(\lambda) = \log L(\hat\theta_\lambda) + 2\,df(\lambda)/n,\qquad
#'       BIC(\lambda) = \log L(\hat\theta_\lambda) + C_n\log(n)\,df(\lambda)/n,}{%
#' AIC = log L + 2 df / n,   BIC = log L + Cn log(n) df / n,}
#' where \eqn{L} is the minimized integrated loss. With `logi = FALSE`,
#' \eqn{AIC = L + 2\,df} and \eqn{BIC = L + C_n\log(n)\,df}.
#' Note that Bottai and Cilluffo (2020) define the log-criteria with
#' \eqn{df/n} and \eqn{\log(n)df/(2n)}; the selected model is the same as with
#' the definition above when \eqn{C_n = 1/2} for the BIC.
#'
#' @return A list with items `call`, `minLambda` (the selected value of
#'   \eqn{\lambda}), `dfMinLambda` (its degrees of freedom), `posMinLambda` (its
#'   position in `object$lambda`), `betaMin` (the coefficients at the selected
#'   \eqn{\lambda}) and `curve` (the criterion along the path).
#' @references
#' Sottile, G., Frumento, P., Chiodi, M. and Bottai, M. (2020). A penalized
#' approach to covariate selection through quantile regression coefficient
#' models. *Statistical Modelling*, 20(4), 369-385. \doi{10.1177/1471082X19825523}
#'
#' Lee, E. R., Noh, H. and Park, B. U. (2014). Model selection via Bayesian
#' information criterion for quantile regression models. *Journal of the
#' American Statistical Association*, 109(505), 216-229.
#' \doi{10.1080/01621459.2013.836975}
#' @seealso [piqr()], [cv_piqr()]
#' @examples
#' set.seed(1)
#' n <- 200
#' x1 <- runif(n); x2 <- runif(n)
#' y <- 1 + qnorm(runif(n)) + 2 * x1
#' fit <- piqr(y ~ x1 + x2, formula.p = ~ slp(p, 2), nlambda = 20)
#' gof.piqr(fit, method = "AIC")
#' gof.piqr(fit, method = "BIC", df = "parameters", plot = FALSE)$posMinLambda
#' @export
gof.piqr <- function(object, method = c("BIC", "AIC"), Cn = 1, plot = TRUE,
                     logi = TRUE, df = NULL, ...) {
  cl <- match.call()
  if ("df.new" %in% names(list(...)))
    stop("the argument 'df.new' was removed in qrcmNP 1.0.0: please use 'df'", call. = FALSE)
  method <- match.arg(method)
  if (is.null(df)) df <- if (identical(object$method, "legacy")) "covariates" else "parameters"
  df.type <- match.arg(df, c("covariates", "parameters"))
  n <- nrow(object$internal$X)
  q <- ncol(object$internal$X)
  if (is.character(Cn)) Cn <- eval(str2lang(Cn), list(n = n, q = q, p = q), baseenv())
  if (!is.numeric(Cn) || length(Cn) != 1 || Cn <= 0) stop("'Cn' must be a positive number", call. = FALSE)
  dfs <- if (df.type == "covariates") vapply(object$seqS, function(s) sum(rowSums(s) > 0), numeric(1))
         else object$npar
  L <- object$obj.function
  curve <- if (logi) {
    switch(method, AIC = log(L) + 2 * dfs / n, BIC = log(L) + Cn * log(n) * dfs / n)
  } else {
    switch(method, AIC = L + 2 * dfs, BIC = L + Cn * log(n) * dfs)
  }
  pos <- which.min(curve)
  if (plot) {
    color <- gray.colors(2)
    plot(log(object$lambda), curve, type = "l", col = color[1], xlab = expression(log(lambda)),
         ylab = method, ...)
    abline(v = log(object$lambda[pos]), col = color[2], lty = 2)
    axis(3, at = log(object$lambda[pos]), labels = paste0("df=", dfs[pos]))
  }
  out <- list(call = cl, minLambda = object$lambda[pos], dfMinLambda = dfs[pos],
              posMinLambda = pos, betaMin = object$coefficients[[pos]], curve = curve,
              method = method)
  class(out) <- "gof.piqr"
  if (plot) invisible(out) else out
}

#' @export
print.gof.piqr <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nSelection of lambda by", x$method, "\n\n")
  cat("Selected lambda:", format(x$minLambda, digits = digits), " (position", x$posMinLambda, ")\n")
  cat("Degrees of freedom:", x$dfMinLambda, "\n\n")
  cat("Coefficients:\n")
  coe <- x$betaMin; coe[coe == 0] <- NA
  print.default(format(coe, digits = digits), print.gap = 2L, quote = FALSE)
  cat("\n")
  invisible(x)
}

#' Unpenalized Re-estimation of the Model Selected by a Penalized QRCM
#'
#' Following Sottile et al. (2020, Section 3.2), the penalized estimator is used
#' to select the covariates, and the selected model is re-estimated without
#' penalty with [qrcm::iqr()], which provides asymptotic standard errors,
#' Wald tests and goodness-of-fit diagnostics.
#'
#' @param object an object of class `"piqr"`.
#' @param pos.lambda the position of the selected value of \eqn{\lambda}. If
#'   missing, it is selected with [gof.piqr()] using `criterion`.
#' @param criterion `"BIC"` (default) or `"AIC"`, used if `pos.lambda` is missing.
#' @param keep `"covariates"` (default): the selected covariates are
#'   re-estimated with all the basis functions allowed by the original `s`;
#'   `"coefficients"`: only the nonzero coefficients are re-estimated.
#' @param data an optional data frame; by default, the `data` of the original call.
#' @param ... further arguments passed to [gof.piqr()].
#' @return An object of class `"iqr"` (see [qrcm::iqr()]), with an additional
#'   attribute `"pos.lambda"`.
#' @note The standard errors are conditional on the selected model and do not
#'   account for the selection step.
#' @seealso [piqr()], [gof.piqr()], [qrcm::iqr()], [qrcm::summary.iqr()]
#' @examples
#' set.seed(1)
#' n <- 300
#' x1 <- runif(n); x2 <- runif(n)
#' dat <- data.frame(y = 1 + qnorm(runif(n)) + 2 * x1, x1, x2)
#' fit <- piqr(y ~ x1 + x2, formula.p = ~ slp(p, 2), data = dat, nlambda = 20)
#' m <- refit_piqr(fit)
#' summary(m)
#' @export
refit_piqr <- function(object, pos.lambda, criterion = c("BIC", "AIC"),
                       keep = c("covariates", "coefficients"), data, ...) {
  criterion <- match.arg(criterion)
  keep <- match.arg(keep)
  if (missing(pos.lambda)) pos.lambda <- gof.piqr(object, method = criterion, plot = FALSE, ...)$posMinLambda
  l <- lambda_index(object, pos.lambda)
  s0 <- object$internal$s
  sel <- object$seqS[[l]]
  snew <- if (keep == "covariates") s0 * (rowSums(sel) > 0 | object$internal$pf == 0)
          else pmax(sel, s0 * (object$internal$pf == 0))
  storage.mode(snew) <- "double"
  dimnames(snew) <- NULL
  cl <- object$call
  keepargs <- intersect(names(cl), c("", "formula", "formula.p", "weights", "data"))
  cl <- cl[keepargs]
  cl[[1L]] <- quote(qrcm::iqr)
  cl$s <- snew
  if (!missing(data)) cl$data <- data
  fit <- eval(cl, parent.frame())
  attr(fit, "pos.lambda") <- l
  fit
}

#' Cross-Validation for Penalized QRCM
#'
#' K-fold cross-validation of the integrated loss along the path of a
#' [piqr()] fit (uncensored data only).
#'
#' @param formula,formula.p,weights,data,s,penalty,penalty.factor,nlambda,lambda.min.ratio,lambda
#'   as in [piqr()].
#' @param nfolds number of folds (default 10).
#' @param foldid an optional vector of fold indicators (integers in
#'   `1:nfolds`), one for each observation used in the fit (i.e., each row of
#'   `data` without missing values in the model variables).
#' @param ... further arguments passed to [piqr()] (e.g., `tol`, `maxit`).
#' @return An object of class `"cv.piqr"`, a list with items `lambda`, `cvm`
#'   (mean cross-validated loss), `cvsd` (its standard error), `lambda.min`,
#'   `lambda.1se`, `pos.min`, `pos.1se` (positions in `lambda`), `foldid` and
#'   `fit` (the [piqr()] fit on the full data).
#' @details The cross-validated loss of observation \eqn{i} is
#'   \eqn{\int_0^1 \rho_p\{y_i - x_i^\top\hat\theta_{-i}\, b(p)\}\,dp},
#'   computed numerically, where \eqn{\hat\theta_{-i}} is estimated without
#'   the fold containing \eqn{i}. In each fold the variables are standardized
#'   as in the full data and \eqn{\lambda} is rescaled by the size of the
#'   training sample, so that the same value of \eqn{\lambda} corresponds to
#'   the same penalty in all fits. `lambda.1se` is the largest \eqn{\lambda}
#'   whose error is within one standard error of the minimum.
#' @seealso [piqr()], [gof.piqr()]
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 300
#' dat <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
#' dat$y <- 1 + qnorm(runif(n)) + 2 * dat$x1
#' cv <- cv_piqr(y ~ x1 + x2 + x3, formula.p = ~ slp(p, 2), data = dat,
#'               nfolds = 5, nlambda = 20)
#' cv
#' plot(cv)
#' coef(cv$fit, pos.lambda = cv$pos.1se)
#' }
#' @export
cv_piqr <- function(formula, formula.p = ~ slp(p, 3), weights, data, s,
                    penalty = c("lasso", "group"), penalty.factor, nlambda = 100,
                    lambda.min.ratio, lambda, nfolds = 10, foldid, ...) {
  penalty <- match.arg(penalty)
  cl <- match.call()
  cl$nfolds <- cl$foldid <- NULL
  cl$penalty <- penalty
  cl[[1L]] <- quote(qrcmNP::piqr)
  full <- eval(cl, parent.frame())
  if (full$internal$type != "iqr") stop("cross-validation is only available for uncensored data", call. = FALSE)
  mf <- full$internal$mf
  n <- nrow(mf)
  lam <- full$lambda
  if (missing(foldid)) foldid <- sample(rep_len(seq_len(nfolds), n))
  if (length(foldid) != n)
    stop(sprintf("'foldid' must have one value for each of the %d observations used in the fit", n), call. = FALSE)
  nfolds <- max(foldid)
  w <- model.weights(mf); if (is.null(w)) w <- rep.int(1, n)

  S <- attr(mf, "stats")
  stats.full <- list(X = S$X, y = S$y, B = S$B)
  mt <- attr(mf, "terms")
  Xall <- model.matrix(mt, mf, contrasts.arg = attr(mf, "contrasts"))
  yall <- model.response(mf)
  pg <- (1:1023) / 1024
  Bp <- apply_bfun(attr(mf, "bfun"), pg, "bfun")
  loss_new <- function(theta, X, y) {
    R <- y - X %*% theta %*% t(Bp)
    rowSums(R * (matrix(pg, nrow(R), ncol(R), byrow = TRUE) - (R < 0))) / 1024
  }
  dots <- list(...)
  tol <- if (is.null(dots$tol)) 1e-6 else dots$tol
  maxit <- if (is.null(dots$maxit)) 100 else dots$maxit
  err <- matrix(NA_real_, n, length(lam))
  for (f in seq_len(nfolds)) {
    test <- which(foldid == f)
    mtr <- mf[-test, , drop = FALSE]
    attr(mtr, "terms") <- mt
    A <- check_in(mtr, full$internal$formula.p, full$internal$s, stats = stats.full)
    ntr <- length(A$V$y)
    fit <- suppressWarnings(piqr_path(A, penalty = penalty, pf = full$internal$pf, nlambda = length(lam),
                                      lambda.min.ratio = NULL, lambda = lam * ntr / n,
                                      tol = tol, maxit = maxit, verbose = FALSE))
    for (j in seq_along(lam)) err[test, j] <- loss_new(fit$theta[[j]], Xall[test, , drop = FALSE], yall[test])
  }
  cvm <- colSums(err * w) / sum(w)
  cvsd <- vapply(seq_along(lam), function(j) {
    fm <- tapply(err[, j] * w, foldid, sum) / tapply(w, foldid, sum)
    stats::sd(fm) / sqrt(nfolds)
  }, numeric(1))
  pos.min <- which.min(cvm)
  pos.1se <- min(which(cvm <= cvm[pos.min] + cvsd[pos.min]))
  out <- list(lambda = lam, cvm = cvm, cvsd = cvsd, lambda.min = lam[pos.min],
              lambda.1se = lam[pos.1se], pos.min = pos.min, pos.1se = pos.1se,
              foldid = foldid, fit = full)
  class(out) <- "cv.piqr"
  out
}

#' @export
print.cv.piqr <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCross-validated integrated loss (", max(x$foldid), " folds)\n\n", sep = "")
  tab <- cbind(Lambda = x$lambda[c(x$pos.min, x$pos.1se)], Index = c(x$pos.min, x$pos.1se),
               Measure = x$cvm[c(x$pos.min, x$pos.1se)], SE = x$cvsd[c(x$pos.min, x$pos.1se)],
               Df = x$fit$df[c(x$pos.min, x$pos.1se)])
  rownames(tab) <- c("min", "1se")
  print(signif(tab, digits))
  cat("\n")
  invisible(x)
}

#' @export
plot.cv.piqr <- function(x, ...) {
  L <- list(...)
  ll <- log(x$lambda)
  args <- list(x = ll, y = x$cvm, ylim = range(x$cvm - x$cvsd, x$cvm + x$cvsd), pch = 20,
               col = "red", xlab = expression(log(lambda)), ylab = "Cross-validated loss")
  args[names(L)] <- L
  do.call(plot, args)
  graphics::segments(ll, x$cvm - x$cvsd, ll, x$cvm + x$cvsd, col = "darkgray")
  abline(v = ll[c(x$pos.min, x$pos.1se)], lty = 3)
  axis(3, at = pretty(ll), labels = approx(ll, x$fit$df, xout = pretty(ll), rule = 2,
                                           method = "constant")$y, tcl = NA)
  invisible(NULL)
}
