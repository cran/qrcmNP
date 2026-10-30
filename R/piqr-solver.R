# Penalized integrated loss minimisation (PILM): proximal Newton solver.
#
# At the current value theta, the integrated loss L is replaced by its
# quadratic approximation  g'(z - theta) + (z - theta)' H (z - theta) / 2,
# where g and H are the exact gradient and Hessian of L (closed form under the
# linear QRCM). The penalised quadratic problem is solved by coordinate descent
# (lasso, soft-thresholding) or block coordinate descent (group lasso, group
# soft-thresholding), and the step z - theta is accepted after a backtracking
# line search on the penalised objective (Lee, Sun and Saunders, 2014). With
# censored or truncated data the estimating equations of Frumento and Bottai
# (2017) are not the gradient of a loss function: the line search then uses the
# norm of the KKT residual as merit function.
#
# All computations are on the internal (standardised) scale: x and b(p) are
# standardised as in Sottile et al. (2020), Section 3, and y is mapped to [0, 10].

#' @noRd
soft <- function(x, a) sign(x) * pmax(abs(x) - a, 0)

#' Value of the penalty.
#' pen: list(type, W (q x k element weights), gw (row weights), free (q x k logical))
#' @noRd
pen_value <- function(theta, pen) {
  if (pen$type == "lasso") sum(pen$W * abs(theta))
  else sum(pen$gw * sqrt(rowSums((theta * pen$free)^2)))
}

#' KKT residual (minimum-norm subgradient of the penalised loss) on the
#' elements in 'idx' (vector positions). Returns a vector.
#' @noRd
kkt_residual <- function(g, theta, lambda, pen, idx) {
  th <- c(theta); r <- numeric(length(g))
  if (pen$type == "lasso") {
    W <- c(pen$W)
    nz <- th != 0
    r[nz] <- g[nz] + lambda * W[nz] * sign(th[nz])
    r[!nz] <- soft(g[!nz], lambda * W[!nz])
  } else {
    q <- nrow(theta); k <- ncol(theta)
    G <- matrix(g, q, k) * pen$free
    TH <- theta * pen$free
    R <- matrix(0, q, k)
    for (j in seq_len(q)) {
      if (!any(pen$free[j, ])) next
      nt <- sqrt(sum(TH[j, ]^2))
      if (pen$gw[j] == 0) R[j, ] <- G[j, ]
      else if (nt > 0) R[j, ] <- G[j, ] + lambda * pen$gw[j] * TH[j, ] / nt
      else {
        ng <- sqrt(sum(G[j, ]^2))
        R[j, ] <- if (ng > 0) G[j, ] * max(0, 1 - lambda * pen$gw[j] / ng) else 0
      }
    }
    r <- c(R)
  }
  r[idx]
}

#' Penalized quadratic subproblem
#'   min_z  gg'(z - th) + (z - th)' H (z - th) / 2 + lambda * pen(z),
#' with pen(z) = sum(W * |z|) (lasso) or sum_g gw_g * ||z_g|| (group lasso;
#' 'groups' is a list of positions). Zero weights = unpenalized coordinates.
#' Coordinate descent identifies the set of nonzero coordinates; the solution
#' is then "polished" exactly on that set (a linear system for the lasso, a few
#' Newton steps for the group lasso) and accepted if it satisfies the KKT
#' conditions. This avoids the slow convergence of coordinate descent when H
#' is ill-conditioned (strongly correlated columns of x (x) b(p)).
#' @noRd
cd_quadratic <- function(gg, H, th, lambda, type, W = NULL, groups = NULL, gw = NULL,
                         tol = 1e-9, maxit = 200) {
  d <- diag(H)
  r0 <- gg - c(H %*% th)                 # gradient of the quadratic at z = 0
  grad_q <- function(z) r0 + c(H %*% z)  # gg + H (z - th)
  if (lambda == 0 || (type == "lasso" && all(W == 0)) || (type == "group" && all(gw == 0))) {
    return(c(solve(H, -r0)))
  }
  if (type == "lasso") {
    kkt_ok <- function(z) {
      gq <- grad_q(z)
      nz <- z != 0
      all(abs(gq[nz] + lambda * W[nz] * sign(z[nz])) <= 1e-7 * (1 + abs(gq[nz]))) &&
        all(abs(gq[!nz]) <= lambda * W[!nz] * (1 + 1e-7) + 1e-10)
    }
    polish <- function(z) {
      A <- which(z != 0 | W == 0)
      if (length(A) == 0) return(z)
      for (rep in 1:5) {
        sg <- sign(z[A]); sg[W[A] == 0] <- 0
        zA <- solve(H[A, A, drop = FALSE], -(r0[A] + lambda * W[A] * sg))
        bad <- W[A] > 0 & sign(zA) != sg
        if (!any(bad)) { z[] <- 0; z[A] <- zA; return(z) }
        # a coefficient changed sign: drop it and retry
        z[A[bad]] <- 0
        A <- A[!bad]
        if (length(A) == 0) { z[] <- 0; return(z) }
      }
      z
    }
    z <- th
    u <- c(H %*% z) + r0                 # gradient of the quadratic at z
    full <- TRUE
    for (it in seq_len(maxit)) {
      set <- if (full) seq_along(z) else which(z != 0 | W == 0)
      maxch <- 0
      for (m in set) {
        zn <- soft(d[m] * z[m] - u[m], lambda * W[m]) / d[m]
        dz <- zn - z[m]
        if (dz != 0) {
          u <- u + H[, m] * dz
          z[m] <- zn
          maxch <- max(maxch, abs(dz) * sqrt(d[m]))
        }
      }
      if (full && it > 1) {
        zp <- polish(z)
        if (kkt_ok(zp)) return(zp)
      }
      if (maxch < tol) { if (full) break; full <- TRUE } else full <- !full && it %% 5 == 0
    }
    zp <- polish(z)
    return(if (kkt_ok(zp)) zp else z)
  }

  # group lasso
  Lg <- vapply(groups, function(G) {
    if (length(G) == 1) H[G, G] else max(eigen(H[G, G], symmetric = TRUE, only.values = TRUE)$values)
  }, numeric(1))
  kkt_ok <- function(z) {
    gq <- grad_q(z)
    for (gi in seq_along(groups)) {
      G <- groups[[gi]]; nz <- sqrt(sum(z[G]^2))
      if (gw[gi] == 0 || nz > 0) {
        r <- gq[G] + (if (gw[gi] > 0) lambda * gw[gi] * z[G] / nz else 0)
        if (max(abs(r)) > 1e-7 * (1 + max(abs(gq[G])))) return(FALSE)
      } else if (sqrt(sum(gq[G]^2)) > lambda * gw[gi] * (1 + 1e-7) + 1e-10) return(FALSE)
    }
    TRUE
  }
  polish <- function(z) {
    act <- which(vapply(seq_along(groups), function(gi) gw[gi] == 0 || any(z[groups[[gi]]] != 0), TRUE))
    if (length(act) == 0) return(z)
    A <- unlist(groups[act])
    for (it in 1:30) {
      gq <- grad_q(z)[A]
      Hn <- H[A, A, drop = FALSE]
      pos <- 0
      for (gi in act) {
        G <- groups[[gi]]; ii <- pos + seq_along(G); pos <- pos + length(G)
        if (gw[gi] > 0) {
          zg <- z[G]; nz <- sqrt(sum(zg^2))
          if (nz == 0) return(z)
          gq[ii] <- gq[ii] + lambda * gw[gi] * zg / nz
          Hn[ii, ii] <- Hn[ii, ii] + lambda * gw[gi] * (diag(length(G)) / nz - tcrossprod(zg) / nz^3)
        }
      }
      step <- tryCatch(solve(Hn, gq), error = function(e) NULL)
      if (is.null(step)) return(z)
      znew <- z; znew[A] <- z[A] - step
      # a group collapsing through zero: stop polishing
      for (gi in act) if (gw[gi] > 0 && sum(znew[groups[[gi]]] * z[groups[[gi]]]) <= 0) return(z)
      z <- znew
      if (max(abs(step)) < tol) break
    }
    z
  }
  z <- th
  u <- c(H %*% z) + r0
  full <- TRUE
  for (it in seq_len(maxit)) {
    set <- if (full) seq_along(groups)
           else which(vapply(groups, function(G) any(z[G] != 0), TRUE) | gw == 0)
    maxch <- 0
    for (gi in set) {
      G <- groups[[gi]]
      zG <- z[G]
      if (gw[gi] == 0) {
        zn <- zG - solve(H[G, G, drop = FALSE], u[G])
      } else {
        v <- zG - u[G] / Lg[gi]
        nv <- sqrt(sum(v^2))
        thr <- lambda * gw[gi] / Lg[gi]
        zn <- if (nv <= thr) 0 * v else (1 - thr / nv) * v
      }
      dz <- zn - zG
      if (any(dz != 0)) {
        u <- u + c(H[, G, drop = FALSE] %*% dz)
        z[G] <- zn
        maxch <- max(maxch, max(abs(dz) * sqrt(d[G])))
      }
    }
    if (full && it > 1) {
      zp <- polish(z)
      if (kkt_ok(zp)) return(zp)
    }
    if (maxch < tol) { if (full) break; full <- TRUE } else full <- !full && it %% 5 == 0
  }
  zp <- polish(z)
  if (kkt_ok(zp)) zp else z
}

#' Make a symmetric positive definite metric from the Jacobian.
#' @noRd
spd_metric <- function(J, ridge = 1e-4) {
  H <- (J + t(J)) / 2
  H <- H + diag(ridge * max(1, mean(abs(diag(H)))), nrow(H))
  for (i in 1:20) {
    if (!inherits(try(chol(H), silent = TRUE), "try-error")) return(H)
    H <- H + diag(10^(i - 4) * max(1, mean(abs(diag(H)))), nrow(H))
  }
  diag(pmax(abs(diag(J)), 1), nrow(J))
}

#' Proximal Newton for a fixed lambda over the active coordinates 'active'
#' (logical q x k). Returns list(theta, ps, g, converged, n.it).
#' @noRd
pnewton <- function(theta, V, bfun, type, lambda, pen, active, tol = 1e-6, maxit = 100,
                    cd.tol = 1e-9, H0 = NULL) {
  ee <- get_ee(type)
  idx <- which(c(active))
  theta[!active] <- 0
  ps <- pstar(theta, V, bfun, type)
  use.loss <- (type == "iqr")
  merit <- function(th, ps, g) {
    if (use.loss) loss_internal(th, V, bfun, type, ps) + lambda * pen_value(th, pen)
    else sum(kkt_residual(g, th, lambda, pen, idx)^2)
  }
  G <- ee(theta, V$y, V$z, V$d, V$X, V$Xw, bfun, ps$y, ps$z, J = TRUE)
  f0 <- merit(theta, ps, G$g)
  converged <- FALSE
  nw <- sum(V$weights)
  best <- list(theta = theta, ps = ps, G = G, f = if (use.loss) Inf else f0)
  it <- 0
  if (length(idx) == 0) return(list(theta = theta, ps = ps, g = G$g, converged = TRUE, n.it = 0))
  if (pen$type == "lasso") {
    cdW <- c(pen$W)[idx]; cdG <- cdgw <- NULL
  } else {
    rows <- ((idx - 1) %% nrow(theta)) + 1
    cdG <- split(seq_along(idx), rows)
    cdgw <- pen$gw[as.integer(names(cdG))]
    cdW <- NULL
  }
  for (it in seq_len(maxit)) {
    g <- G$g
    Jx <- G$J[idx, idx, drop = FALSE]
    # at theta = 0 (models without unpenalized coefficients) all the p* are at the
    # extremes of the grid and J = 0: use the metric H0 instead
    if (!is.null(H0) && max(abs(Jx)) < 1e-8) Jx <- H0[idx, idx, drop = FALSE]
    H <- spd_metric(Jx)
    z <- cd_quadratic(g[idx], H, theta[idx], lambda, pen$type, W = cdW, groups = cdG, gw = cdgw,
                      tol = cd.tol)
    delta <- z - theta[idx]
    if (max(abs(delta)) < tol) { converged <- TRUE; break }
    th.new <- theta
    pv0 <- pen_value(theta, pen)
    tstep <- 1
    accepted <- FALSE
    for (ls in 1:40) {
      th.new[idx] <- theta[idx] + tstep * delta
      ps.new <- pstar(th.new, V, bfun, type)
      if (use.loss) {
        f.new <- merit(th.new, ps.new, NULL)
        Delta <- sum(g[idx] * delta) + lambda * (pen_value({tmp <- theta; tmp[idx] <- z; tmp}, pen) - pv0)
        if (f.new <= f0 + 1e-4 * tstep * min(Delta, 0)) { accepted <- TRUE; break }
      } else {
        G.new <- ee(th.new, V$y, V$z, V$d, V$X, V$Xw, bfun, ps.new$y, ps.new$z, J = FALSE)
        f.new <- merit(th.new, ps.new, G.new$g)
        if (f.new < f0) { accepted <- TRUE; break }
      }
      tstep <- tstep / 2
    }
    if (!accepted) {
      if (use.loss) {
        # no decrease is possible at the precision of the internal grid
        converged <- (max(abs(delta)) < sqrt(tol))
        break
      }
      # estimating equations: take the full step (fixed-point iteration)
      tstep <- 1
      th.new[idx] <- theta[idx] + delta
      ps.new <- pstar(th.new, V, bfun, type)
      G.new <- ee(th.new, V$y, V$z, V$d, V$X, V$Xw, bfun, ps.new$y, ps.new$z, J = FALSE)
      f.new <- merit(th.new, ps.new, G.new$g)
    }
    theta <- th.new; ps <- ps.new; f0 <- f.new
    G <- ee(theta, V$y, V$z, V$d, V$X, V$Xw, bfun, ps$y, ps$z, J = TRUE)
    if (!use.loss && f0 < best$f) best <- list(theta = theta, ps = ps, G = G, f = f0)
    if (max(abs(tstep * delta)) < tol) { converged <- TRUE; break }
    if (!use.loss && sqrt(f0) / nw < 1e-5) { converged <- TRUE; break }
  }
  if (!use.loss) {
    # estimating equations: return the iterate with the smallest KKT residual;
    # with censored data the equations are piecewise constant (internal grid of
    # p), and an exact zero cannot be attained.
    theta <- best$theta; ps <- best$ps; G <- best$G
    kres <- max(abs(kkt_residual(G$g, theta, lambda, pen, idx))) / nw
    if (getOption("qrcmNP.debug", FALSE)) message("KKT residual / n: ", signif(kres, 3))
    converged <- converged || (kres < 5e-4)
  }
  list(theta = theta, ps = ps, g = G$g, converged = converged, n.it = it)
}
