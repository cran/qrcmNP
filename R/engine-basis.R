# Internal engine: basis functions b(p), bisection, numerical integration.
#
# Most functions in this file are adapted from the package 'qrcm'
# (Paolo Frumento, GPL-2, https://CRAN.R-project.org/package=qrcm), so that
# 'qrcmNP' does not need to access non-exported objects of 'qrcm'.

#' @noRd
pmax0 <- function(x) (x + abs(x)) / 2

#' Row-wise Kronecker (Khatri-Rao) product used to vectorise theta-related
#' quantities. Column (h - 1) * q + j of the result is the product of column j
#' of X and column h of B, which
#' matches the column-major vectorisation of the q x k matrix theta.
#' @noRd
kr <- function(X, B) {
  q <- ncol(X); k <- ncol(B)
  X[, rep.int(seq_len(q), k), drop = FALSE] * B[, rep(seq_len(k), each = q), drop = FALSE]
}

#' Evaluate b(p) ("bfun") or b'(p) ("b1fun") for the external basis.
#' @noRd
apply_bfun <- function(bfun, p, fun = c("bfun", "b1fun")) {
  fun <- match.arg(fun)
  k <- attr(bfun, "k")
  n <- length(p)
  if (inherits(bfun, "slp.basis")) {
    pp <- matrix(NA_real_, n, k + 1)
    pp[, 1] <- 1; pp[, 2] <- p
    if (k > 1) for (j in 2:k) pp[, j + 1] <- pp[, j] * p
    out <- if (fun == "bfun") tcrossprod(pp, t(bfun$a))
           else cbind(0, tcrossprod(pp[, 1:k, drop = FALSE], t(bfun$a1[-1, -1, drop = FALSE])))
    if (!attr(bfun, "intB")) out <- out[, -1, drop = FALSE]
  } else {
    out <- matrix(NA_real_, n, k)
    if (fun == "bfun") for (j in 1:k) out[, j] <- bfun[[j]](p)
    else for (j in 1:k) out[, j] <- bfun[[j]](p, deriv = 1)
  }
  out
}

#' Bisection for external use (precision ~ 1e-6). Returns F(y | x).
#' @noRd
p_bisec <- function(theta, y, X, bfun, n.it = 20) {
  n <- length(y); k <- ncol(theta)
  bp <- attr(bfun, "bp")
  p <- attr(bfun, "p")
  Xtheta <- tcrossprod(X, t(theta))
  eta <- tcrossprod(Xtheta, bp[512, , drop = FALSE])
  m <- as.integer(512 + sign(y - eta) * 256)
  for (i in 3:10) {
    eta <- .rowSums(Xtheta * bp[m, , drop = FALSE], n, k)
    m <- m + as.integer(sign(y - eta) * (2^(10 - i)))
  }
  m <- p[m]
  for (i in 11:n.it) {
    bp <- apply_bfun(bfun, m, "bfun")
    delta.m <- y - .rowSums(Xtheta * bp, n, k)
    m <- m + sign(delta.m) / 2^i
  }
  m <- c(m)
  out.l <- which(m == 1 / 2^n.it)
  out.r <- which(m == 1 - 1 / 2^n.it)
  m[out.l] <- 0; m[out.r] <- 1
  attr(m, "out") <- c(out.l, out.r)
  attr(m, "out.r") <- out.r
  m
}

#' Bisection on the index scale of the internal 1023-point grid (precision
#' ~ 0.001). Returns integer indices with attributes "out" and "out.r".
#' @noRd
p_bisec_internal <- function(theta, y, X, bp) {
  n <- length(y); k <- ncol(theta)
  Xtheta <- tcrossprod(X, t(theta))
  eta <- tcrossprod(Xtheta, bp[512, , drop = FALSE])
  m <- as.integer(512 + sign(y - eta) * 256)
  for (i in 3:10) {
    eta <- .rowSums(Xtheta * bp[m, , drop = FALSE], n, k)
    m <- m + as.integer(sign(y - eta) * (2^(10 - i)))
  }
  attr(m, "out") <- c(which(m == 1), which(m == 1023))
  attr(m, "out.r") <- which(m == 1023)
  m
}

#' Monotone (if possible) interpolating spline of a basis function.
#' @noRd
make_bfun <- function(p, x) {
  n <- length(x)
  x1 <- x[1:(n - 1)]; x2 <- x[2:n]
  method <- if (all(x1 < x2) || all(x1 > x2)) "hyman" else "fmm"
  stats::splinefun(p, x, method = method)
}

#' Cumulative trapezoidal integral ("int") or forward difference ("der") of
#' the columns of fx on a grid with increments dx.
#' @noRd
num_fun <- function(dx, fx, op = c("int", "der")) {
  op <- match.arg(op)
  n <- length(dx) + 1
  fL <- fx[1:(n - 1), , drop = FALSE]
  fR <- fx[2:n, , drop = FALSE]
  if (op == "int") {
    out <- rbind(0, 0.5 * dx * (fL + fR))
    out <- apply(out, 2, cumsum)
    if (!is.matrix(out)) out <- matrix(out, nrow = n)
  } else {
    out <- (fR - fL) / dx
    out <- rbind(out[1, ], out)
  }
  out
}

#' Shifted Legendre polynomial basis (coefficients of the polynomial, of its
#' derivative and of its first and second integrals).
#' @noRd
slp_basis <- function(k, intercept) {
  K <- k + 1
  a <- matrix(0, K, K)
  for (i1 in 0:k) for (i2 in 0:i1) a[i2 + 1, i1 + 1] <- choose(i1, i2) * choose(i1 + i2, i2)
  a[, seq.int(2, K, 2)] <- -a[, seq.int(2, K, 2)]
  a[seq.int(2, K, 2), ] <- -a[seq.int(2, K, 2), ]
  a1 <- A <- AA <- matrix(NA_real_, K, K)
  for (j in 0:k) {
    a1[j + 1, ] <- a[j + 1, ] * j
    A[j + 1, ] <- a[j + 1, ] / (j + 1)
    AA[j + 1, ] <- A[j + 1, ] / (j + 2)
  }
  if (!intercept) a[1, -1] <- A[1, -1] <- AA[1, -1] <- 0
  out <- list(a = a, a1 = a1, A = A, AA = AA)
  attr(out, "k") <- k
  class(out) <- "slp.basis"
  out
}

#' Is formula.p exactly ~ slp(p, k) (possibly without intercept)?
#' @noRd
is_slp <- function(f) {
  test.p <- seq(0, 1, 0.1)
  B <- stats::model.matrix(f, data = data.frame(p = test.p))
  if (nrow(B) == 0) return(FALSE)
  a <- attr(B, "assign")
  if (any(a > 1)) return(FALSE)
  B <- B[, a == 1, drop = FALSE]
  k <- ncol(B)
  intercept <- FALSE
  if (any(B != qrcm::slp(test.p, k = k, intercept))) {
    intercept <- TRUE
    if (any(B != qrcm::slp(test.p, k = k, intercept))) return(FALSE)
  }
  out <- TRUE
  attr(out, "k") <- k
  attr(out, "intercept") <- intercept
  attr(out, "intB") <- any(a == 0)
  out
}

#' beta(p) = theta b(p) with delta-method standard errors at a single p.
#' model: list(coefficients, covar, mf) where attr(mf, "bfun") is the basis.
#' @noRd
extract_p <- function(model, p, cov = FALSE) {
  theta <- model$coefficients
  v <- model$covar
  q <- nrow(theta); k <- ncol(theta)
  bfun <- attr(model$mf, "bfun")
  pred.p <- apply_bfun(bfun, p, "bfun")
  beta <- c(pred.p %*% t(theta))
  if (is.null(v)) {
    out <- cbind(Estimate = beta)
    rownames(out) <- rownames(theta)
    return(list(coef = out))
  }
  cov.beta <- matrix(NA_real_, q, q)
  for (j1 in 1:q) {
    w1 <- seq.int(j1, q * k, q)
    for (j2 in j1:q) {
      w2 <- seq.int(j2, q * k, q)
      cov.beta[j1, j2] <- cov.beta[j2, j1] <- pred.p %*% v[w1, w2, drop = FALSE] %*% t(pred.p)
    }
  }
  se <- sqrt(diag(cov.beta))
  z <- beta / se
  out <- cbind(beta, se, z, 2 * stats::pnorm(-abs(z)))
  colnames(out) <- c("Estimate", "std.err", "z value", "p(>|z|)")
  rownames(out) <- colnames(cov.beta) <- rownames(cov.beta) <- rownames(theta)
  if (cov) list(coef = out, cov = cov.beta) else list(coef = out)
}

#' beta(p) on a grid of p, as a named list of data frames.
#' @noRd
pred_beta <- function(model, p, se = FALSE) {
  theta <- model$coefficients
  if (se && !is.null(model$covar)) {
    Beta <- SE <- matrix(NA_real_, length(p), nrow(theta))
    for (i in seq_along(p)) {
      b <- extract_p(model, p[i])$coef
      Beta[i, ] <- b[, 1]; SE[i, ] <- b[, 2]
    }
    out <- lapply(seq_len(nrow(theta)), function(j)
      data.frame(p = p, beta = Beta[, j], se = SE[, j],
                 low = Beta[, j] - 1.96 * SE[, j], up = Beta[, j] + 1.96 * SE[, j]))
  } else {
    beta <- apply_bfun(attr(model$mf, "bfun"), p, "bfun") %*% t(theta)
    out <- lapply(seq_len(nrow(theta)), function(j) data.frame(p = p, beta = beta[, j]))
  }
  names(out) <- rownames(theta)
  out
}
