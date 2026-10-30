# The algorithm of 'qrcmNP' 0.2.x (method = "legacy"), ported to the internal
# engine of the package. It is kept for reproducibility of the results in
# Sottile et al. (2020). A covariate (row of theta) enters the model when the
# l1-norm of its gradient exceeds lambda; within the active set, the
# estimating equations  g(theta) = sign * lambda * w  are solved by
# Newton-Raphson, where 'sign' is the sign of the gradient when the covariate
# enters and w = |g_jh| / sum_h |g_jh|. Covariates never leave the model.

#' Newton-Raphson for the (sign-fixed) penalized estimating equations.
#' Port of piqr.newton() of qrcmNP 0.2.2 (itself adapted from qrcm::iqr.newton).
#' @noRd
legacy_newton <- function(theta, V, bfun, s, type, tol, maxit, safeit, eps0, segno = 1, lambda = 0) {
  ee0 <- get_ee(type)
  pl <- c(segno * lambda)
  ee <- function(th, ps, J = FALSE, G) {
    if (missing(G)) {
      out <- ee0(th, V$y, V$z, V$d, V$X, V$Xw, bfun, ps$y, ps$z, J = J)
      out$g <- out$g - pl
      out
    } else ee0(th, V$y, V$z, V$d, V$X, V$Xw, bfun, ps$y, ps$z, J = J, G = G)
  }
  q <- nrow(theta); k <- ncol(theta)
  s <- c(s == 1)
  ps <- pstar(theta, V, bfun, type)
  G <- ee(theta, ps)
  g <- G$g[s]
  conv <- FALSE
  eps <- eps0

  for (i in seq_len(safeit)) {
    if (conv || max(abs(g)) < tol) break
    u <- rep.int(0, q * k); u[s] <- g
    delta <- matrix(u, q, k); delta[is.na(delta)] <- 0
    cond <- FALSE
    while (!cond) {
      new.theta <- theta - delta * eps
      if (max(abs(delta * eps)) < tol) { conv <- TRUE; break }
      ps <- pstar(new.theta, V, bfun, type)
      G1 <- ee(new.theta, ps)
      g1 <- G1$g[s]
      cond <- (sum(g1^2) < sum(g^2))
      eps <- eps * 0.5
    }
    if (conv) break
    g <- g1; G <- G1; theta <- new.theta
    eps <- min(eps * 2, 0.1)
  }

  alg <- "nr"
  conv <- FALSE
  eps <- eps * 10
  ps <- pstar(theta, V, bfun, type)
  h <- ee(theta, ps, J = TRUE, G = G)$J[s, s, drop = FALSE]
  h <- h + diag(0.0001, nrow(h))
  i <- 0
  for (i in seq_len(maxit)) {
    if (conv || max(abs(g)) < tol) break
    if (type == "iqr") {
      H1 <- try(chol(h), silent = TRUE); err <- inherits(H1, "try-error")
    } else {
      H1 <- qr(h); err <- (H1$rank != ncol(h))
    }
    if (!err) {
      if (alg == "gs") { alg <- "nr"; eps <- 1 }
      delta <- if (type == "iqr") chol2inv(H1) %*% g else qr.solve(H1) %*% g
    } else {
      if (alg == "nr") { alg <- "gs"; eps <- 1 }
      delta <- g
    }
    u <- rep.int(0, q * k); u[s] <- delta
    delta <- matrix(u, q, k); delta[is.na(delta)] <- 0
    cond <- FALSE
    while (!cond) {
      new.theta <- theta - delta * eps
      if (max(abs(delta * eps)) < tol) { conv <- TRUE; break }
      ps <- pstar(new.theta, V, bfun, type)
      G1 <- ee(new.theta, ps)
      g1 <- G1$g[s]
      cond <- (sum(g1^2) < sum(g^2))
      eps <- eps * 0.5
    }
    if (conv) break
    g <- g1; G <- G1; theta <- new.theta
    h <- ee(theta, ps, J = TRUE, G = G)$J[s, s, drop = FALSE]
    eps <- if (i > 1) min(eps * 10, 1) else min(eps * 10, 0.1)
  }
  ps <- pstar(theta, V, bfun, type)
  grad <- ee0(theta, V$y, V$z, V$d, V$X, V$Xw, bfun, ps$y, ps$z, J = FALSE)$g
  list(coefficients = matrix(theta, q, k), converged = (i < maxit), n.it = i, ps = ps,
       ee = g, jacobian = h, rank = (alg == "nr") * sum(s), grad = grad)
}

#' Port of divide.et.impera.2() of qrcmNP 0.2.2.
#' @noRd
legacy_divide <- function(fit, V, bfun, s, type, tol, maxit, safeit, eps0, segno, lambda) {
  maxind <- function(A) {
    w <- which.max(A); n <- nrow(A)
    row <- w %% n; if (row == 0) row <- n
    c(row, floor(w / n) + (row < n))
  }
  if (sum(s) == 1) return(fit)
  theta0 <- theta <- fit$coefficients
  ee <- s; ee[s == 1] <- abs(fit$ee)
  failed <- TRUE
  while (failed) {
    i <- maxind(ee)
    sbis <- s * (ee > 0); sbis[-i[1], ] <- 0
    theta <- legacy_newton(theta, V, bfun, sbis, type, tol, 2 * sum(sbis), 5, eps0, segno, lambda)$coefficients
    sbis <- s * (ee > 0); sbis[, -i[2]] <- 0
    theta <- legacy_newton(theta, V, bfun, sbis, type, tol, 2 * sum(sbis), 5, eps0, segno, lambda)$coefficients
    failed <- all(theta == theta0)
    ee[i[1], i[2]] <- 0
    if (all(ee == 0)) break
  }
  if (failed) return(fit)
  legacy_newton(theta, V, bfun, s, type, tol, maxit, safeit, eps0, segno, lambda)
}

#' Port of piqr.internal() of qrcmNP 0.2.2: fit for a given active set.
#' @noRd
legacy_fit <- function(A, s, tol, maxit, segno = 1, lambda = 0) {
  V <- A$V; bfun <- A$internal.bfun; type <- A$type
  maxit <- max(10 + 10 * sum(s), maxit)
  theta0 <- start_theta(V, bfun, s, type)
  safeit <- 10; try.count <- 0; eeTol <- 1; eps0 <- 0.05
  fit.ok <- FALSE
  while (!fit.ok) {
    try.count <- try.count + 1
    fit <- legacy_newton(theta0, V, bfun, s, type, tol, maxit, safeit, eps0, segno, lambda)
    if (max(abs(fit$ee)) > eeTol && try.count > 2)
      fit <- legacy_divide(fit, V, bfun, s, type, tol, maxit, safeit, eps0, segno, lambda)
    fit.ok <- (fit$rank == ncol(fit$jacobian) && max(abs(fit$ee)) < eeTol)
    if (fit.ok || try.count >= 20) break
    eeTol <- eeTol * 2; safeit <- safeit + 2; eps0 <- eps0 / 2
  }
  if (!fit.ok) stop("unable to fit the model", call. = FALSE)
  S <- list(B = A$stats.B, X = A$stats.X, y = A$stats.y)
  obj <- loss_internal(fit$coefficients, V, bfun, type, fit$ps) * (S$y$M - S$y$m) / 10
  covar <- cov_theta(fit$coefficients, V, bfun, fit$ps, type, s)
  out <- check_out(fit$coefficients, S, covar = covar$Q)
  list(theta = out$theta, theta.int = fit$coefficients, covar = out$covar, obj = obj,
       grad = fit$grad, converged = fit$converged)
}

#' Port of the path algorithm of piqr() in qrcmNP 0.2.2.
#' @noRd
legacy_path_core <- function(A, nlambda, lambda.min.ratio, lambda, tol, maxit, verbose) {
  s <- A$s; q <- nrow(s); k <- ncol(s)
  intercept <- A$stats.X$intercept
  cX <- A$stats.X$const
  V <- A$V; bfun <- A$internal.bfun; type <- A$type
  thetas <- thetas.int <- covars <- list(); obj <- dl <- conv <- NULL

  if (!intercept) {
    th0 <- matrix(0, q, k)
    ps <- pstar(th0, V, bfun, type)
    grad <- get_ee(type)(th0, V$y, V$z, V$d, V$X, V$Xw, bfun, ps$y, ps$z, J = FALSE)$g
    s_int <- NULL
  } else {
    s_int <- s[cX, ]
    s1 <- matrix(0, q, k); s1[cX, ] <- s_int
    f <- legacy_fit(A, s1, tol, maxit)
    thetas[[1]] <- f$theta; thetas.int[[1]] <- f$theta.int; covars[[1]] <- f$covar; obj <- f$obj; conv <- f$converged
    grad <- f$grad
  }
  grad[c(s) == 0] <- 0
  dl <- cbind(grad)
  GRAD <- matrix(grad, q, k)
  segno <- sign(GRAD)
  temp <- rowSums(abs(GRAD))
  percent <- abs(GRAD) / temp
  if (is.null(lambda)) {
    l <- max(temp); il <- which.max(temp)
    lmin <- lambda.min.ratio * l
    if (nlambda > 1) l <- exp(seq(log(l), log(lmin), length.out = nlambda))
  } else {
    l <- lambda; nlambda <- length(l); il <- which(temp >= l[1])
  }
  AA <- matrix(FALSE, q, k)
  for (i in il) AA[i, s[i, ] == 1] <- TRUE
  if (intercept) AA[cX, ] <- (s_int == 1)

  if (nlambda == 1 || !intercept) { thetas <- thetas.int <- covars <- list(); obj <- conv <- NULL; dl <- NULL }
  np <- if (nlambda == 1 || !intercept) 0 else 1
  flag <- TRUE
  for (kk in 1:10000) {
    if (flag) { np <- np + 1; if (np > nlambda) { np <- np - 1; break } }
    else {
      thetas <- thetas[-np]; thetas.int <- thetas.int[-np]; covars <- covars[-np]; obj <- obj[-np]; conv <- conv[-np]; dl <- dl[, -np, drop = FALSE]
    }
    flag <- TRUE
    al <- l[np]
    s.np <- 1 * AA
    tempSegno <- s.np * segno
    if (intercept) tempSegno[cX, ] <- 0
    f <- try(legacy_fit(A, s.np, tol, maxit, segno = tempSegno, lambda = al * percent), silent = TRUE)
    if (inherits(f, "try-error")) { np <- np - 1; break }
    thetas[[np]] <- f$theta; thetas.int[[np]] <- f$theta.int; covars[[np]] <- f$covar; obj[np] <- f$obj; conv[np] <- f$converged
    grad <- f$grad; grad[c(s) == 0] <- 0
    dl <- cbind(dl, grad)
    GRAD <- matrix(grad, q, k)
    segno2 <- sign(GRAD)
    temp <- rowSums(abs(GRAD))
    percent <- abs(GRAD) / temp
    if (nlambda != 1) {
      aa <- AA[, 1]
      for (m in seq_len(q)) {
        if (!aa[m] && abs(temp[m]) >= al) {
          flag <- FALSE
          AA[m, ] <- (s[m, ] == 1)
          segno[m, ] <- segno2[m, ]
        }
      }
    } else if (kk > 1) break
    if (flag && verbose)
      message(sprintf("lambda %3d: %10.5f   obj.function = %10.4f   df = %d", np, al, obj[np],
                      sum(rowSums(AA) > 0)))
    if (nlambda == 1 && flag) break
  }
  dl <- apply(dl, 2, function(x) rowSums(abs(matrix(x, q, k))))
  list(lambda = l[seq_len(np)], theta = thetas, theta.int = thetas.int, covar = covars, obj.function = obj,
       dl = matrix(dl, nrow = q), converged = conv)
}

#' @noRd
piqr_legacy <- function(A, nlambda, lambda.min.ratio, lambda, tol, maxit, verbose) {
  legacy_path_core(A, nlambda, lambda.min.ratio, lambda, tol, maxit, verbose)
}
