# Internal engine: data preparation, scaling and back-transformation.
#
# Adapted from 'qrcm' (Paolo Frumento, GPL-2): check.in.iqr() and check.out().
# Differences with respect to 'qrcm':
#  * no principal-component rotation of x and b(p): the penalty of piqr() must
#    act on (standardised) coefficients that keep their meaning, as in
#    Sottile et al. (2020), Section 3;
#  * interval-censored data are not supported.

#' Prepare data, basis and scaling for estimation. If 'stats' is supplied
#' (list(X, y, B) as returned in stats.X, stats.y, stats.B), its location-scale
#' values are used instead of those computed from the data.
#' @return list(mf, U, V, stats.B, stats.X, stats.y, internal.bfun, bfun, s, type)
#' @noRd
check_in <- function(mf, formula.p, s, stats = NULL) {

  if (!any(all.vars(formula.p) == "p"))
    stop("the quantile function must depend on 'p'", call. = FALSE)
  if (!missing(s) && !is.null(s) && all(s == 0))
    stop("'s' cannot be all zero", call. = FALSE)

  explore_s <- function(s, dim) {
    if (dim == 2) s <- t(s)
    out <- 1
    if ((r <- nrow(s)) > 1) {
      for (j in 2:r) {
        done <- FALSE; rj <- s[j, ]
        for (h in 1:(j - 1)) if (all(rj == s[h, ])) { out[j] <- out[h]; done <- TRUE }
        if (!done) out[j] <- max(out) + 1
      }
    }
    out
  }

  # weights
  weights <- stats::model.weights(mf)
  if (!is.null(weights) && any(weights < 0)) stop("negative 'weights'", call. = FALSE)
  if (is.null(weights)) {
    weights <- rep.int(1, nrow(mf)); alarm <- FALSE
  } else {
    alarm <- (weights == 0)
    sel <- which(!alarm)
    mf <- mf[sel, , drop = FALSE]
    weights <- weights[sel]
    weights <- weights / mean(weights)
  }
  if (any(alarm)) warning("observations with null weight will be dropped", call. = FALSE)
  if ((n <- nrow(mf)) == 0) stop("zero non-NA cases", call. = FALSE)

  # y, z, d
  zyd <- stats::model.response(mf)
  type <- attributes(zyd)$type
  zyd <- cbind(zyd)
  if (is.null(type)) {
    y <- zyd[, 1]; z <- rep.int(-Inf, n); d <- rep.int(1, n); type <- fittype <- "iqr"
  } else if (type == "right") {
    y <- zyd[, 1]; z <- rep.int(-Inf, n); d <- zyd[, 2]
    type <- if (any(d == 0)) "ciqr" else "iqr"
    fittype <- "ciqr"
  } else if (type == "counting") {
    z <- zyd[, 1]; y <- zyd[, 2]; d <- zyd[, 3]; type <- fittype <- "ctiqr"
    if (all(z < min(y))) type <- if (any(d == 0)) "ciqr" else "iqr"
  } else {
    stop("only uncensored, 'right' and 'counting' Surv() data are supported", call. = FALSE)
  }
  if (!any(d == 1)) stop("all observations are censored", call. = FALSE)
  attr(type, "fittype") <- fittype

  # x and b(p)
  X <- stats::model.matrix(attr(mf, "terms"), mf); q <- ncol(X)
  termlabelsX <- attr(attr(mf, "terms"), "term.labels")
  assignX <- attr(X, "assign")
  coefnamesX <- colnames(X)

  # p1: grid to build the splinefuns (more points on the tails);
  # p2: grid for external bisection (p reachable by bisection on the p scale);
  # p3: grid for internal use (reachable by bisection on the index scale).
  p1 <- stats::pbeta(seq.int(stats::qbeta(1e-6, 2, 2), stats::qbeta(1 - 1e-6, 2, 2), length.out = 1000), 2, 2)
  p2 <- (1:1023) / 1024
  p3 <- stats::pbeta(seq.int(stats::qbeta(1 / (1000 * n), 2.5, 2.5),
                             stats::qbeta(1 - 1 / (1000 * n), 2.5, 2.5), length.out = 1023), 2.5, 2.5)

  if ((use.slp <- is_slp(formula.p))) {
    k <- attr(use.slp, "k")
    intercept <- attr(use.slp, "intercept")
    intB <- attr(use.slp, "intB")
    assignB <- (1 - intB):k
    termlabelsB <- paste0("slp", 1:k)
    coefnamesB <- if (intB) c("(Intercept)", termlabelsB) else termlabelsB
    k <- k + intB
  } else {
    B <- stats::model.matrix(formula.p, data = data.frame(p = c(p1, p2, p3)))
    B1 <- B[1:1000, , drop = FALSE]
    B2 <- B[1001:2023, , drop = FALSE]
    B3 <- B[2024:3046, , drop = FALSE]
    k <- ncol(B)
    assignB <- attr(B, "assign")
    termlabelsB <- attr(stats::terms(formula.p), "term.labels")
    coefnamesB <- colnames(B)
  }
  if (missing(s) || is.null(s)) {
    s <- matrix(1, q, k)
  } else {
    s <- as.matrix(s)
    if (any(dim(s) != c(q, k))) stop(sprintf("wrong size of 's': it must be a %d x %d matrix", q, k), call. = FALSE)
    if (any(s != 0 & s != 1)) stop("'s' can only contain 0 and 1", call. = FALSE)
  }

  # x singularities (set s = 0 where they occur)
  vx <- qr(X); selx <- vx$pivot[seq_len(vx$rank)]
  if (vx$rank < q) s[-selx, ] <- 0

  # b(p) singularities, dropped row by row based on s
  if (!use.slp && qr(B3)$rank < k) {
    u <- explore_s(s, 1)
    for (j in unique(u)) {
      sel <- which(s[which(u == j)[1], ] == 1)
      if (length(sel) > 1) {
        vbj <- qr(B3[, sel, drop = FALSE])
        if ((rj <- vbj$rank) < length(sel)) s[u == j, sel[-vbj$pivot[1:rj]]] <- 0
      }
    }
  }

  # location-scale statistics for x, b(p), y
  ry <- range(y[is.finite(y)]); my <- ry[1]; My <- ry[2]
  sX <- apply(X, 2, stats::sd); mX <- colMeans(X)
  intX <- (length((constX <- which(sX == 0 & mX != 0))) > 0)
  varsX <- which(sX > 0); zeroX <- which(sX == 0 & mX == 0)
  sX[constX] <- X[1, constX]; mX[constX] <- 0; sX[zeroX] <- 1
  if (length(constX) > 1) { zeroX <- c(zeroX, constX[-1]); constX <- constX[1] }

  if (!use.slp) {
    sB <- apply(B3, 2, stats::sd); mB <- colMeans(B3)
    intB <- (length((constB <- which(sB == 0 & mB != 0))) > 0); varsB <- which(sB > 0)
    if (length(varsB) == 0) stop("the quantile function must depend on 'p'", call. = FALSE)
    if (length(constB) > 1) stop("remove multiple constant functions from 'formula.p'", call. = FALSE)
    if (any(sB == 0 & mB == 0)) stop("remove zero functions from 'formula.p'", call. = FALSE)
    sB[constB] <- B3[1, constB]; mB[constB] <- 0
  } else {
    sB <- rep.int(1, k); mB <- rep.int(0, k)
    if (intB) { constB <- 1; varsB <- 2:k } else { constB <- integer(0); varsB <- 1:k }
  }

  if (all(s[, varsB] == 0))
    stop("the quantile function must depend on 'p' (wrong specification of 's')", call. = FALSE)
  if (!((intX & intB) && s[constX, constB] == 1)) { my <- 0; My <- stats::sd(y[is.finite(y)]) * 5 }
  if (!intX || all(s[constX, ] == 0)) mX <- rep.int(0, q)
  else for (j in varsX) if (any(s[j, ] > s[constX, ])) mX[j] <- 0
  if (!intB || any(s[, constB] == 0)) mB <- rep.int(0, k)
  if (!is.null(stats)) {
    # reuse the scaling of another fit (e.g., the full data in cross-validation),
    # so that the same lambda corresponds to the same penalty
    mX <- stats$X$m; sX <- stats$X$s; my <- stats$y$m; My <- stats$y$M
    if (!use.slp) { mB <- stats$B$m; sB <- stats$B$s }
  }

  # external bfun (post-estimation)
  if (!use.slp) {
    bfun <- list()
    if (intB) bfun[[constB]] <- function(p, deriv = 0) rep.int(1 - deriv, length(p))
    for (j in varsB) bfun[[j]] <- make_bfun(p1, B1[, j])
    names(bfun) <- coefnamesB
    attr(bfun, "k") <- k
  } else {
    bfun <- slp_basis(k - intB, intercept)
    if (!intB) bfun$a[1, 1] <- bfun$A[1, 1] <- bfun$AA[1, 1] <- 0
    attr(bfun, "intB") <- intB
    B2 <- apply_bfun(bfun, p2, "bfun")
  }
  attr(bfun, "bp") <- B2
  attr(bfun, "p") <- p2

  # scaling of x, b(p), y (no rotations)
  U <- list(X = X, y = y, z = z)
  X <- scale(X, center = mX, scale = sX)
  attr(X, "scaled:center") <- attr(X, "scaled:scale") <- NULL
  y <- (y - my) / (My - my) * 10
  z <- z0 <- (z - my) / (My - my) * 10
  z[z < min(y)] <- -Inf
  if (!use.slp) {
    B3 <- scale(B3, center = mB, scale = sB)
    attr(B3, "scaled:center") <- attr(B3, "scaled:scale") <- NULL
  }
  rotX <- diag(1, q); MX <- rep.int(0, q); SX <- rep.int(1, q)
  rotB <- diag(1, k); MB <- rep.int(0, k); SB <- rep.int(1, k)

  # pre-evaluated internal basis on grid p3
  p <- p3
  if (!use.slp) {
    bp <- B3
    dp <- p[-1] - p[-1023]
    b1p <- num_fun(dp, bp, "der")
    Bp <- num_fun(dp, bp, "int")
    BBp <- num_fun(dp, Bp, "int")
    BB1 <- BBp[1023, ]
  } else {
    kk <- attr(bfun, "k")
    pp <- matrix(NA_real_, 1023, kk + 1)
    pp[, 1] <- 1; pp[, 2] <- p
    if (kk > 1) for (j in 2:kk) pp[, j + 1] <- pp[, j] * p
    bp <- tcrossprod(pp, t(bfun$a))
    b1p <- cbind(0, tcrossprod(pp[, 1:kk, drop = FALSE], t(bfun$a1[-1, -1, drop = FALSE])))
    pp <- cbind(pp, pp[, kk + 1] * p, pp[, kk + 1] * p^2)
    Bp <- tcrossprod(pp[, 2:(kk + 2)], t(bfun$A))
    BBp <- tcrossprod(pp[, 3:(kk + 3)], t(bfun$AA))
    BB1 <- colSums(bfun$AA)
    if (!intB) {
      bp <- bp[, -1, drop = FALSE]; b1p <- b1p[, -1, drop = FALSE]
      Bp <- Bp[, -1, drop = FALSE]; BBp <- BBp[, -1, drop = FALSE]
      BB1 <- BB1[-1]
    }
  }
  BB1 <- matrix(rep(BB1, each = n), n)
  internal.bfun <- list(p = p, bp = bp, b1p = b1p, Bp = Bp, BBp = BBp, BB1 = BB1)
  attr(internal.bfun, "pfun") <- stats::approxfun(c(p[1], 0.5 * (p[-1023] + p[-1])), p,
                                                  method = "constant", rule = 2)

  stats.B <- list(m = mB, s = sB, M = MB, S = SB, rot = rotB, const = constB, vars = varsB,
                  intercept = intB, term.labels = termlabelsB, assign = assignB, coef.names = coefnamesB)
  stats.X <- list(m = mX, s = sX, M = MX, S = SX, rot = rotX, const = constX, vars = varsX,
                  intercept = intX, term.labels = termlabelsX, assign = assignX, coef.names = coefnamesX)
  stats.y <- list(m = my, M = My)

  V <- list(X = X, Xw = X * weights, y = y, z = z, d = d, weights = weights)
  if (type == "ctiqr") V$z0 <- z0
  dimnames(s) <- list(coefnamesX, coefnamesB)
  list(mf = mf, U = U, V = V, stats.B = stats.B, stats.X = stats.X, stats.y = stats.y,
       internal.bfun = internal.bfun, bfun = bfun, s = s, type = type)
}

#' Back-transform theta (and its covariance) from the internal to the original scale.
#' @noRd
check_out <- function(theta, S, covar = NULL) {

  blockdiag <- function(A, d, type = 1) {
    h <- nrow(A); g <- d / h
    if (type == 1) {
      out <- diag(1, d)
      for (j in 1:g) { ind <- (j * h - h + 1):(j * h); out[ind, ind] <- A }
    } else {
      out <- matrix(0, d, d)
      for (i1 in 1:h) for (i2 in 1:h) {
        ind1 <- (i1 * g - g + 1):(i1 * g)
        ind2 <- (i2 * g - g + 1):(i2 * g)
        out[ind1, ind2] <- diag(A[i1, i2], g)
      }
      out <- t(out)
    }
    out
  }
  mydiag <- function(x) if (length(x) > 1) diag(x) else matrix(x, 1, 1)

  th <- cbind(c(theta))
  q <- nrow(theta); k <- ncol(theta); g <- q * k
  if (is.null(covar)) covar <- matrix(0, g, g)
  aX <- S$X; ay <- S$y; aB <- S$B
  cX <- aX$const; cB <- aB$const

  tr <- function(A) { th <<- A %*% th; covar <<- A %*% covar %*% t(A) }

  tr(blockdiag(mydiag(1 / aX$S), g))
  if (aX$intercept) { A <- diag(1, q); A[cX, ] <- -aX$M; A[cX, cX] <- 1; tr(blockdiag(A, g)) }
  tr(blockdiag(mydiag(1 / aB$S), g, 2))
  if (aB$intercept) { A <- diag(1, k); A[, cB] <- -aB$M; A[cB, cB] <- 1; tr(blockdiag(A, g, 2)) }
  tr(blockdiag(aX$rot, g))
  tr(blockdiag(t(aB$rot), g, 2))
  tr(blockdiag(mydiag(1 / aX$s), g))
  if (aX$intercept) { A <- diag(1, q); A[cX, ] <- -aX$m / aX$s[cX]; A[cX, cX] <- 1; tr(blockdiag(A, g)) }
  tr(blockdiag(mydiag(1 / aB$s), g, 2))
  if (aB$intercept) { A <- diag(1, k); A[, cB] <- -aB$m / aB$s[cB]; A[cB, cB] <- 1; tr(blockdiag(A, g, 2)) }

  v <- (ay$M - ay$m) / 10
  th <- th * v
  covar <- covar * (v^2)
  theta <- matrix(th, q, k)
  theta[cX, cB] <- theta[cX, cB] + ay$m / aB$s[cB] / aX$s[cX]
  list(theta = theta, covar = covar)
}
