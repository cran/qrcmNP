# Internal engine: integrated loss, gradient (estimating equations), Jacobian,
# covariance and starting values for the linear QRCM
#   Q(p | x, theta) = x' theta b(p).
#
# Adapted from 'qrcm' (Paolo Frumento, GPL-2): iqr.ee(), ciqr.ee(), ctiqr.ee(),
# iobjfun(), iobjfun.ct(), cov.fun.iqr(), start.iqr(), km().
# The Jacobians are computed with a single cross-product of row-wise Kronecker
# (Khatri-Rao) matrices instead of k^2 separate cross-products.
#
# Notation: p.star.y are integer indices on the internal grid bfun$p, with
# attributes "out" (at the extremes of the grid) and "out.r" (right extreme).

#' Gradient and Jacobian for uncensored data.
#' g = X' W (BB1 - B(p*)), J = Z' diag(w f) Z with Z = b(p*) (x) x and f = 1/Q'(p*).
#' @noRd
iqr_ee <- function(theta, y, z, d, X, Xw, bfun, p.star.y, p.star.z, J = TRUE, G, i = FALSE) {
  n <- length(y)
  k <- ncol(theta)
  if (missing(G)) {
    B <- bfun$Bp[p.star.y, , drop = FALSE]
    S1 <- bfun$BB1 - B
    g <- if (!i) c(crossprod(Xw, S1)) else kr(X, S1)
  } else {
    B <- G$B; g <- G$g
  }
  if (J) {
    b <- bfun$bp[p.star.y, , drop = FALSE]
    b1 <- bfun$b1p[p.star.y, , drop = FALSE]
    A1 <- 1 / c(.rowSums(tcrossprod(X, t(theta)) * b1, n, k))
    A1 <- pmax0(A1)
    A1[attr(p.star.y, "out")] <- 0
    A1[!is.finite(A1)] <- 0
    J <- crossprod(kr(Xw, b) * A1, kr(X, b))
  }
  list(g = g, J = J, B = B)
}

#' Gradient and Jacobian for right-censored data.
#' @noRd
ciqr_ee <- function(theta, y, z, d, X, Xw, bfun, p.star.y, p.star.z, J = TRUE, G, i = FALSE) {
  n <- length(y)
  k <- ncol(theta)
  BB1 <- bfun$BB1
  if (missing(G)) {
    B <- bfun$Bp[p.star.y, , drop = FALSE]
    BB <- bfun$BBp[p.star.y, , drop = FALSE]
    py <- bfun$p[p.star.y]
    a <- (1 - d) / (1 - py); a[attr(p.star.y, "out.r")] <- 0
    S1 <- a * (BB - BB1)
    a <- d; a[attr(p.star.y, "out.r")] <- 1
    S1 <- BB1 - a * B + S1
    g <- if (!i) c(crossprod(Xw, S1)) else kr(X, S1)
  } else {
    B <- G$B; BB <- G$BB; g <- G$g; py <- G$py
  }
  if (J) {
    b <- bfun$bp[p.star.y, , drop = FALSE]
    b1 <- bfun$b1p[p.star.y, , drop = FALSE]
    A1 <- 1 / c(.rowSums(tcrossprod(X, t(theta)) * b1, n, k))
    A1 <- pmax0(A1)
    A1[attr(p.star.y, "out")] <- 0
    A1[!is.finite(A1)] <- 0
    a <- 1 - py
    a[attr(p.star.y, "out.r")] <- 1 # arbitrary, cancelled by the zeroes in A1
    A2 <- d * b + (1 - d) / (a^2) * (BB1 - B * a - BB)
    J <- crossprod(kr(Xw, A2) * A1, kr(X, b))
  }
  list(g = g, J = J, B = B, BB = BB, py = py)
}

#' Gradient and Jacobian for censored and truncated data.
#' @noRd
ctiqr_ee <- function(theta, y, z, d, X, Xw, bfun, p.star.y, p.star.z, J = TRUE, G, i = FALSE) {
  n <- length(y)
  k <- ncol(theta)
  BB1 <- bfun$BB1
  if (missing(G)) {
    B.y <- bfun$Bp[p.star.y, , drop = FALSE]
    BB.y <- bfun$BBp[p.star.y, , drop = FALSE]
    B.z <- bfun$Bp[p.star.z, , drop = FALSE]
    BB.z <- bfun$BBp[p.star.z, , drop = FALSE]
    py <- bfun$p[p.star.y]
    pz <- bfun$p[p.star.z]
    out.y <- attr(p.star.y, "out.r")
    out.z <- attr(p.star.z, "out.r")
    a <- d
    S1.y <- (1 - d) / (1 - py) * (BB.y - BB1)
    S1.z <- 1 / (1 - pz) * (BB.z - BB1)
    S1.y[out.y, ] <- 0; a[out.y] <- 1
    S1.y[out.z, ] <- S1.z[out.z, ] <- a[out.z] <- 0
    S1 <- -a * B.y + S1.y - S1.z
    g <- if (!i) c(crossprod(Xw, S1)) else kr(X, S1)
  } else {
    B.y <- G$B.y; BB.y <- G$BB.y; B.z <- G$B.z; BB.z <- G$BB.z; g <- G$g; py <- G$py; pz <- G$pz
  }
  if (J) {
    b.y <- bfun$bp[p.star.y, , drop = FALSE]
    b1.y <- bfun$b1p[p.star.y, , drop = FALSE]
    b.z <- bfun$bp[p.star.z, , drop = FALSE]
    b1.z <- bfun$b1p[p.star.z, , drop = FALSE]
    Xtheta <- tcrossprod(X, t(theta))
    A1.y <- pmax0(1 / c(.rowSums(Xtheta * b1.y, n, k)))
    A1.z <- pmax0(1 / c(.rowSums(Xtheta * b1.z, n, k)))
    A1.y[attr(p.star.y, "out")] <- 0
    A1.z[attr(p.star.z, "out")] <- 0
    A1.y[!is.finite(A1.y)] <- 0; A1.z[!is.finite(A1.z)] <- 0
    ay <- 1 - py; az <- 1 - pz
    ay[attr(p.star.y, "out.r")] <- az[attr(p.star.z, "out.r")] <- 1
    A2.y <- d * b.y - (1 - d) / (ay^2) * (B.y * ay + BB.y - BB1)
    A2.z <- 1 / (az^2) * (B.z * az + BB.z - BB1)
    J <- crossprod(kr(Xw, A2.y * A1.y), kr(X, b.y)) + crossprod(kr(Xw, A2.z * A1.z), kr(X, b.z))
  }
  list(g = g, J = J, B.y = B.y, BB.y = BB.y, B.z = B.z, BB.z = BB.z, py = py, pz = pz)
}

#' Select the estimating equation for a given data type.
#' @noRd
get_ee <- function(type) {
  switch(as.character(type), iqr = iqr_ee, ciqr = ciqr_ee, ctiqr = ctiqr_ee,
         stop("unsupported data type", call. = FALSE))
}

#' CDF indices (p*) for y and, with truncation, z.
#' @noRd
pstar <- function(theta, V, bfun, type) {
  out <- list(y = p_bisec_internal(theta, V$y, V$X, bfun$bp), z = NULL)
  if (type == "ctiqr") out$z <- p_bisec_internal(theta, V$z, V$X, bfun$bp)
  out
}

#' Integrated loss, uncensored data (internal scale):
#' sum_i w_i { y_i (p*_i - 1/2) + x_i' theta (BB1 - B(p*_i)) }.
#' @noRd
iobjfun <- function(theta, y, X, weights, bfun, p.star) {
  B <- bfun$Bp[p.star, , drop = FALSE]
  py <- bfun$p[p.star]
  S <- sweep(-B, 2, bfun$BB1[1, ], "+")
  sum(y * weights * (py - 0.5)) + sum(((X * weights) %*% theta) * S)
}

#' Integrated loss for censored and truncated data (on a 250-point subgrid),
#' as in 'qrcm'. Note: this is not the function whose gradient defines the
#' estimating equations of ciqr_ee()/ctiqr_ee().
#' @noRd
iobjfun_ct <- function(theta, z, y, d, X, weights, bfun, py, pz, type) {
  trunc <- (type == "ctiqr")
  n <- nrow(X)
  r <- 250
  psel <- round(seq.int(1, 1023, length.out = r))
  p <- bfun$p[psel]
  bp <- bfun$bp[psel, , drop = FALSE]
  dp <- p[-1] - p[-r]; dp <- c(dp[1], dp)
  beta <- tcrossprod(theta, bp)
  eta <- tcrossprod(X, t(beta))
  deltay <- y - eta
  omegay <- (deltay <= 0)
  Sy <- 1 - matrix(pmin(py, 1 - 1e-6), n, r)
  deltaz <- if (trunc) z - eta else -Inf
  omegaz <- if (trunc) deltaz <= 0 else 1
  Sz <- 1 - (if (trunc) matrix(pmin(pz, 1 - 1e-6), n, r) else 0)
  p <- t(matrix(t(p), r, n))
  dp <- t(matrix(t(dp), r, n))
  pbar <- 1 - p
  omega.hat <- omegay * (1 - (1 - d) * pbar / Sy) - omegaz * (1 - pbar / Sz)
  loss <- -weights * (omega.hat * deltay)
  sum(.rowSums(loss * dp, n, r))
}

#' Loss on the internal scale for any data type.
#' @noRd
loss_internal <- function(theta, V, bfun, type, ps) {
  if (type == "iqr") return(iobjfun(theta, V$y, V$X, V$weights, bfun, ps$y))
  py <- bfun$p[ps$y]
  pz <- if (type == "ctiqr") {
    pz <- pmin(bfun$p[ps$z], py - 1e-8); pz[ps$z == 1] <- 0; pz
  } else NULL
  iobjfun_ct(theta, V$z, V$y, V$d, V$X, V$weights, bfun, py, pz, type)
}

#' Sandwich covariance of theta (internal scale), as in qrcm::cov.fun.iqr.
#' @noRd
cov_theta <- function(theta, V, bfun, ps, type, s) {
  ee <- get_ee(type)
  s <- c(s == 1)
  G.i <- ee(theta, V$y, V$z, V$d, V$X, V$Xw, bfun, ps$y, ps$z, J = TRUE, i = TRUE)
  s.i <- G.i$g[, s, drop = FALSE]
  G <- G.i$J[s, s, drop = FALSE]
  # regularisation inherited from 'qrcm' (>= 3.0), to avoid singularities
  s.i <- t(t(s.i) - colMeans(s.i))
  G <- G + diag(0.01, ncol(G))
  Omega <- chol2inv(chol(crossprod(s.i * V$weights, s.i) + diag(0.01, ncol(s.i))))
  Q0 <- t(G) %*% Omega %*% G + diag(0.01, ncol(G))
  Q <- chol2inv(chol(Q0))
  if ((cc <- ncol(Q)) > 1) {
    sdq <- sqrt(diag(Q))
    bound <- 0.999 * tcrossprod(sdq)
    diag(bound) <- Inf
    Q <- sign(Q) * pmin(abs(Q), bound)
  }
  U <- matrix(0, length(s), length(s))
  U[s, s] <- Q
  list(Q = U, Q0 = Q0)
}

#' Kaplan-Meier (or Lynden-Bell) estimate of the marginal CDF evaluated at y.
#' @noRd
km_cdf <- function(z, y, d, w, type) {
  fit <- tryCatch({
    m <- if (type == "ctiqr") {
      survival::survfit(survival::coxph(survival::Surv(z, y, d) ~ 1, weights = pmax(w, 1e-6),
                                        timefix = FALSE), type = "kaplan-meier")
    } else {
      survival::survfit(survival::Surv(y, d) ~ 1, weights = w)
    }
    list(time = m$time, cdf = 1 - m$surv)
  }, error = function(e) NULL)
  if (is.null(fit)) return((rank(y, ties.method = "first") - 0.5) / length(y))
  out <- stats::approx(fit$time, fit$cdf, xout = y, method = "constant", rule = 2, f = 0)$y
  pmin(pmax(out, 0.5 / length(y)), 1 - 0.5 / length(y))
}

#' Starting values: least squares of y on x (x) b(p*), with p* estimated
#' marginally (ranks, or Kaplan-Meier with censoring/truncation), restricted to
#' the free parameters in s. Robust trimming as in qrcm::start.iqr.
#' @noRd
start_theta <- function(V, bfun, s, type) {
  y <- V$y
  p.star <- if (type == "iqr") (rank(y, ties.method = "first") - 0.5) / length(y)
            else km_cdf(V$z, y, V$d, V$weights, type)
  p.star <- attr(bfun, "pfun")(p.star)
  b.star <- bfun$bp[match(p.star, bfun$p), , drop = FALSE]
  Z <- kr(V$X, b.star)[, c(s) == 1, drop = FALSE]
  w <- V$weights
  restol <- 4
  repeat {
    m <- stats::lm.wfit(Z, y, w)
    res <- m$residuals
    ok <- abs(res) / stats::sd(res) < restol
    if (all(ok) || sum(ok) < 0.5 * length(y)) break
    Z <- Z[ok, , drop = FALSE]; y <- y[ok]; w <- w[ok]
    restol <- restol + 1
  }
  out <- rep.int(0, length(s))
  out[c(s) == 1] <- m$coefficients
  out[is.na(out)] <- 0
  matrix(out, nrow(s))
}
