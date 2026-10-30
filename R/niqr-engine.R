# Internal engine of niqr(): integrated loss, gradient and Hessian of a
# nonlinear parametric quantile model (Bottai and Cilluffo, 2020)
#
#   Q(p | x_i, theta) = fun2(fun(theta, p), x_i),
#
# where fun(theta, p) returns the matrix of "coefficients" beta(p | theta)
# (one row for each value of p) and fun2(beta, X) combines them row by row with
# the covariates. If fun2 is NULL the model is linear in the covariates,
# Q = x' beta(p | theta).
#
# The integral over p is computed with the trapezoidal rule on a grid that is
# denser in the tails, splitting the cell that contains p*_i = F(y_i | x_i),
# where the integrand is not differentiable. Loss and gradient are therefore
# consistent with each other and smooth in theta. The Hessian is the
# "expected" (Gauss-Newton) Hessian
#   H = sum_i w_i f_i dQ_i(p*_i) dQ_i(p*_i)',  f_i = 1 / Q'(p*_i | x_i),
# which is positive semi-definite and consistent for the Hessian of the
# expected loss at the true parameter (the term involving second derivatives
# of Q has zero expectation at the true value).

#' Grid on (0, 1) for numerical integration.
#' @noRd
niqr_grid <- function(n.p, eps) {
  p <- stats::pbeta(seq.int(stats::qbeta(eps, 2, 2), stats::qbeta(1 - eps, 2, 2), length.out = n.p), 2, 2)
  h <- diff(p)
  list(p = p, h = h, w = c(h, 0) / 2 + c(0, h) / 2, eps = eps)
}

#' Quantile function and its derivatives at paired (p_i, x_i), or at a single p
#' for all observations. Returns list(Q, D) with D = dQ/dtheta (n x d) if deriv.
#' @noRd
niqr_qd <- function(theta, p, M, deriv = FALSE) {
  X <- M$X; n <- nrow(X)
  ev <- function(th) {
    b <- M$fun(th, p)
    if (!is.matrix(b)) b <- matrix(b, nrow = length(p))
    if (is.null(M$fun2)) {
      if (length(p) == 1) c(X %*% c(b)) else .rowSums(b * X, n, ncol(X))
    } else {
      if (length(p) == 1) b <- b[rep.int(1L, n), , drop = FALSE]
      c(M$fun2(b, X))
    }
  }
  Q <- ev(theta)
  if (!deriv) return(list(Q = Q))
  d <- length(theta)
  D <- matrix(NA_real_, n, d)
  for (k in seq_len(d)) {
    h <- M$h * max(1, abs(theta[k]))
    tp <- tm <- theta; tp[k] <- tp[k] + h; tm[k] <- tm[k] - h
    D[, k] <- (ev(tp) - ev(tm)) / (2 * h)
  }
  list(Q = Q, D = D)
}

#' Quantile function on the whole grid (n x m matrix).
#' @noRd
niqr_Qgrid <- function(theta, M) {
  pg <- M$grid$p
  if (is.null(M$fun2)) {
    b <- M$fun(theta, pg)
    if (!is.matrix(b)) b <- matrix(b, nrow = length(pg))
    return(tcrossprod(M$X, b))
  }
  vapply(pg, function(pp) niqr_qd(theta, pp, M)$Q, numeric(nrow(M$X)))
}

#' p*_i = F(y_i | x_i, theta) by bracketing on the grid and bisection.
#' Returns p* with attribute "cell" (0 = below the grid, m = above the grid).
#' @noRd
niqr_pstar <- function(theta, M, Qg) {
  pg <- M$grid$p; m <- length(pg); y <- M$y
  cell <- .rowSums(Qg <= y, nrow(Qg), m)
  lo <- pg[pmax(cell, 1)]; hi <- pg[pmin(cell + 1, m)]
  inside <- which(cell > 0 & cell < m)
  ps <- numeric(length(y))
  ps[cell == 0] <- 0; ps[cell == m] <- 1
  if (length(inside)) {
    a <- lo[inside]; b <- hi[inside]
    Mi <- M; Mi$X <- M$X[inside, , drop = FALSE]
    yi <- y[inside]
    for (it in seq_len(M$n.bisec)) {
      mid <- (a + b) / 2
      up <- niqr_qd(theta, mid, Mi)$Q <= yi
      a[up] <- mid[up]; b[!up] <- mid[!up]
    }
    ps[inside] <- (a + b) / 2
  }
  attr(ps, "cell") <- cell
  ps
}

#' Integrated loss (per observation, before weighting).
#' @noRd
niqr_loss_i <- function(theta, M, Qg = niqr_Qgrid(theta, M), ps = niqr_pstar(theta, M, Qg)) {
  g <- M$grid; pg <- g$p; m <- length(pg)
  R <- M$y - Qg
  Fm <- (matrix(pg, nrow(R), m, byrow = TRUE) - (R < 0)) * R      # check loss at the nodes (>= 0)
  L <- c(Fm %*% g$w)
  # exact split of the cell containing p*, where the integrand has a kink
  ins <- split_cells(ps, R, m)
  if (length(ins)) {
    c0 <- attr(ps, "cell")[ins]
    f0 <- Fm[cbind(ins, c0)]; f1 <- Fm[cbind(ins, c0 + 1)]
    L[ins] <- L[ins] - (f0 + f1) * g$h[c0] / 2 +
      f0 * (ps[ins] - pg[c0]) / 2 + f1 * (pg[c0 + 1] - ps[ins]) / 2
  }
  attr(L, "ps") <- ps
  L
}

#' Observations whose residual changes sign once in the cell containing p*
#' (always the case if Q is increasing in p).
#' @noRd
split_cells <- function(ps, R, m) {
  cell <- attr(ps, "cell")
  ins <- which(cell > 0 & cell < m)
  if (!length(ins)) return(ins)
  ok <- R[cbind(ins, cell[ins])] >= 0 & R[cbind(ins, cell[ins] + 1)] < 0
  ins[ok]
}

#' Loss, individual gradients, gradient and Gauss-Newton Hessian.
#' @noRd
niqr_eval <- function(theta, M, hessian = TRUE) {
  g <- M$grid; pg <- g$p; m <- length(pg)
  n <- nrow(M$X); d <- length(theta)
  Qg <- niqr_Qgrid(theta, M)
  ps <- niqr_pstar(theta, M, Qg)
  L <- niqr_loss_i(theta, M, Qg, ps)
  cell <- attr(ps, "cell")
  R <- M$y - Qg
  inside <- seq_len(n) %in% split_cells(ps, R, m)

  # gradient: - int (p - I(p >= p*)) dQ/dtheta dp
  S <- matrix(0, n, d)
  D0 <- D1 <- matrix(0, n, d)
  lin <- is.null(M$fun2)
  if (lin) {
    # dbeta/dtheta on the grid: m x q x d
    q <- ncol(M$X)
    Jb <- array(NA_real_, c(m, q, d))
    for (k in seq_len(d)) {
      h <- M$h * max(1, abs(theta[k]))
      tp <- tm <- theta; tp[k] <- tp[k] + h; tm[k] <- tm[k] - h
      bp <- M$fun(tp, pg); bm <- M$fun(tm, pg)
      Jb[, , k] <- matrix(bp - bm, nrow = m) / (2 * h)
    }
  }
  for (l in seq_len(m)) {
    Dl <- if (lin) M$X %*% matrix(Jb[l, , ], q, d) else niqr_qd(theta, pg[l], M, deriv = TRUE)$D
    coefl <- pg[l] - (R[, l] < 0)
    S <- S - g$w[l] * coefl * Dl
    i0 <- which(cell == l); i1 <- which(cell == l - 1)
    if (length(i0)) D0[i0, ] <- Dl[i0, , drop = FALSE]
    if (length(i1)) D1[i1, ] <- Dl[i1, , drop = FALSE]
  }
  # derivatives at p* (observations inside the grid only)
  ii <- which(inside)
  Ds <- matrix(0, n, d)
  if (length(ii)) {
    Mi <- M; Mi$X <- M$X[ii, , drop = FALSE]
    Ds[ii, ] <- niqr_qd(theta, ps[ii], Mi, deriv = TRUE)$D
    c0 <- cell[ii]; hc <- g$h[c0]
    p0 <- pg[c0]; p1 <- pg[c0 + 1]; pst <- ps[ii]
    f0 <- -p0 * D0[ii, , drop = FALSE]                       # node c0, left of p*: I = 0
    f1 <- -(p1 - 1) * D1[ii, , drop = FALSE]                 # node c0 + 1, right of p*: I = 1
    fl <- -pst * Ds[ii, , drop = FALSE]                      # at p*-
    fr <- -(pst - 1) * Ds[ii, , drop = FALSE]                # at p*+
    S[ii, ] <- S[ii, ] - (f0 + f1) * hc / 2 + (f0 + fl) * (pst - p0) / 2 + (fr + f1) * (p1 - pst) / 2
  }
  w <- M$weights
  out <- list(loss = sum(w * L) / sum(w), L.i = L, S = S, grad = colSums(w * S) / sum(w),
              ps = ps, cell = cell)
  if (hessian) {
    f <- pdf <- rep.int(0, n)
    if (length(ii)) {
      dl <- 1e-5
      pa <- pmin(ps[ii] + dl, 1 - 1e-12); pb <- pmax(ps[ii] - dl, 1e-12)
      Qp <- (niqr_qd(theta, pa, Mi)$Q - niqr_qd(theta, pb, Mi)$Q) / (pa - pb)
      pdf[ii] <- 1 / Qp
      f[ii] <- pmax0(1 / Qp)
      f[!is.finite(f)] <- 0
    }
    out$H <- crossprod(Ds * (w * f), Ds) / sum(w)
    out$pdf <- pdf
  }
  out
}
