#' @importFrom stats splinefun sd printCoefmat runif pnorm predict
#' @importFrom graphics plot points abline par axis
#' @importFrom grDevices n2mfrow
#' @importFrom utils menu setTxtProgressBar txtProgressBar

Q <- function(beta, X){
  # beta <- fun(theta, p)
  X <- as.matrix(X)
  q <- ncol(X)
  out <- 0
  for(i in 1:q) out <- out + beta[, i]*X[, i]
  out
}

funQ1 <- function(X){
  X <- as.matrix(X)
  return(X)
}

boot <- function(fun, fun2, X, y, start_theta, true_theta, p, seed, boot=FALSE, R=100, ...){
  if(!missing(seed)) set.seed(seed)
  if(missing(p)) p <- seq(.01, .99, .01)
  if(missing(start_theta)) stop("Please insert start values!")
  if(missing(true_theta)) stop("Please insert true theta values!")

  Xb <- X
  yb <- y

  n <- nrow(X)
  q <- ncol(X)
  np <- length(p)
  ntheta <- length(true_theta)
  true_beta <- fun(true_theta, p)
  control <- list(...)
  if(missing(fun2)) fun2 <- Q

  lowTheta <- upTheta <- seTheta <- theta <- matrix(NA, nrow=ntheta, ncol=R, dimnames=list(paste0("theta", 0:(ntheta-1)), 1:R))
  lowBeta <- upBeta <- seBeta <- beta <- array(NA, dim=c(R, np, q), dimnames=list(1:R, p, c("(Intercept)", paste0("X", 1:(q-1)))))

  pb <- txtProgressBar(min=0, max=R, style=3)
  for(i in 1:R){
    setTxtProgressBar(pb, i)

    if(boot){
      id <- sample(n, n, replace=T)
      Xb <- Xb[id, ]
      yb <- yb[id]
    }else{
      bet <- fun(true_theta, runif(n))
      yb <- fun2(bet, X=Xb)
    }

    o <- try(niqr(fun, fun2, start_theta, X=Xb, y=yb, control=control), silent=TRUE)
    if(class(o) != "try-error"){
      theta[, i] <- o$x
      seTheta[, i] <- o$se
      lowTheta[, i] <- o$x - 1.96*o$se
      upTheta[, i] <- o$x + 1.96*o$se

      o2 <- predict(o, type="beta", p=p)
      for(j in 1:q){
        beta[i, , j] <- o2[[j]][, 2]
        seBeta[i, , j] <- o2[[j]][, 3]
        lowBeta[i, , j] <- o2[[j]][, 4]
        upBeta[i, , j] <- o2[[j]][, 5]
      }
    }
  }
  close(pb)

  matTheta <- cbind(true_theta=true_theta, fit_theta=apply(theta, 1, mean, na.rm=T))
  matSeTheta <- cbind(fit_se=rowMeans(seTheta, na.rm=T), boot_se=apply(theta, 1, sd, na.rm=T))
  copSeTheta <- NULL
  for(i in 1:R) copSeTheta <- cbind(copSeTheta, ((true_theta >= lowTheta[, i]) & (true_theta <= upTheta[, i])))

  matBeta <- cbind(true_beta, apply(beta, 3, function(.x) apply(.x, 2, mean, na.rm=TRUE)))
  matSeBeta <- cbind(apply(seBeta, 3, function(.x) colMeans(.x, na.rm=TRUE)), apply(beta, 3, function(.x) apply(.x, 2, sd, na.rm=TRUE)))
  copSeBeta <- array(NA, dim=c(R, np, q), dimnames=list(1:R, p, c("(Intercept)", paste0("X", 1:(q-1)))))
  for(i in 1:R){
    for(j in 1:(q)){
      copSeBeta[i, , j] <- ((true_beta[, j] >= lowBeta[i, , j]) & (true_beta[, j] <= upBeta[i, , j]))
    }
  }


  return(list(theta=theta, seTheta=seTheta, lowTheta=lowTheta, upTheta=upTheta, matTheta=matTheta, matSeTheta=matSeTheta, copSeTheta=copSeTheta,
              beta=beta, seBeta=seBeta, lowBeta=lowBeta, upBeta=upBeta, matBeta=matBeta, matSeBeta=matSeBeta, copSeBeta=copSeBeta))
}

check.fun <- function(obj, b){
  if(!missing(b)) obj$x <- b
  n.var <- length(obj$x)
  n <- n2mfrow(n.var)
  par(mfrow=n)
  for(j in 1:n.var){
    seqX <- seq(obj$x[j]-1, obj$x[j]+1, by=.01)
    seqY <- sapply(seqX, function(i, j){
      bb <- obj$x
      bb[j] <- i
      Fun(obj$internal$fun, obj$internal$fun2, bb, obj$internal$y, obj$internal$X, q=obj$internal$q, n=obj$internal$n, p0=obj$internal$p0, npp=obj$internal$npp)}, j=j)
    plot(seqX, seqY, type="l")
    abline(v=obj$x[j])
    if(missing(b)) abline(v=c(obj$x[j]-1.96*obj$se[j], obj$x[j]+1.96*obj$se[j]), col=3)
    abline(v=seqX[which.min(seqY)], col=2)
    axis(3, at=seqX[which.min(seqY)], labels=round(seqX[which.min(seqY)],2))
  }
  par(mfrow=c(1,1))
}

#' @export
print.niqr <- function (x, digits = max(3L, getOption("digits") - 3L), ...){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  q <- ncol(x$internal$X)
  s <- x$s
  s2 <- c(0, cumsum(s))
  seCoe <- coe <- matrix(NA, q, length(x$x), dimnames=list(c("(Intercept)", if(q > 1) paste0("X", 1:(q-1))), paste0("theta", 0:(length(x$x)-1))))
  for(j in 1:q){
    coe[j, (s2[j]+1):s2[j+1]] <- x$x[(s2[j]+1):s2[j+1]]
    seCoe[j, (s2[j]+1):s2[j+1]] <- x$se[(s2[j]+1):s2[j+1]]
  }

  cat("Coefficients:\n")
  print.default(format(coe, digits = digits), print.gap = 2L, quote = FALSE)
  cat("\n")

  cat("Standard Errors:\n")
  print.default(format(seCoe, digits = digits), print.gap = 2L, quote = FALSE)
  cat("\n")

  invisible(x)
}

#' @export
summary.niqr <- function(object, p, ...){
  Beta <- list()
  fun <- object$internal$fun
  if(!missing(p)){
    q <- object$internal$q
    bet <- object$internal$fun(object$x, p)
    gr <- grad2(fun, x=object$x, p=p, n.var=object$internal$n.var, q=q, s=object$s, ind=object$ind, h=object$internal$h1, method=object$internal$meth)$grad
    if(length(p) == 1) gr <- t(gr)
    s2 <- c(0, cumsum(object$s))
    coe <- matrix(0, nrow=q, ncol=4, dimnames=list(colnames(object$internal$X), c("Estimate", "std.err", "z value", "p(>|z|)")))
    for(j in 1:length(p)){
      coe[, 1] <- bet[j, ]
      for(i in 1:q){
        coe[i, 2] <- sqrt(t(gr[j, (s2[i]+1):s2[i+1]]) %*% object$covar[(s2[i]+1):s2[i+1], (s2[i]+1):s2[i+1]] %*% as.matrix(gr[j, (s2[i]+1):s2[i+1]]))
      }
      coe[, 3] <- coe[, 1]/coe[, 2]
      coe[, 4] <- 2*(1-pnorm(abs(coe[, 3])))

      Beta[[j]] <- coe
    }
    Beta$p <- p
  }else{
    qq <- length(object$x)
    coe <- matrix(0, nrow=qq, ncol=4, dimnames=list(paste0("theta", 0:(qq-1)), c("Estimate", "std.err", "z value", "p(>|z|)")))
    coe[, 1] <- object$x
    coe[, 2] <- object$se
    coe[, 3] <- coe[, 1]/coe[, 2]
    coe[, 4] <- 2*(1-pnorm(abs(coe[,3])))

    Beta$coe <- coe
  }

  Beta$obj <- object
  class(Beta) <- "summary.niqr"
  return(Beta)
}

#' @export
print.summary.niqr <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  cat("\nCall: ", paste(deparse(x$obj$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if(length(x) <= 2){

    cat("n. of iterations:", x$obj$iter, "\n")
    cat("n. of observations:", x$obj$internal$n, "\n\n")

    cat("######################", "\n")
    cat("######################", "\n\n")

    cat("Coefficients:\n")
    printCoefmat(x$coe, digits=max(3L, getOption("digits") - 3L), signif.stars=TRUE, cs.ind=1:2, tst.ind=3, P.values=TRUE, has.Pvalue=TRUE)
    cat("\n")

    cat("######################", "\n")
    cat("######################", "\n")

    if(!is.null(x$obj$fx)){
      cat("\n")
      cat("Minimized loss function:", x$obj$fx)
    }
  }

  else{
    for(j in 1:length(x$p)){
      cat("\n", paste("p = ", x$p[j], "\n"))
      cat("\n")
      cat("Coefficients:\n")
      printCoefmat(x[[j]], digits=digits, signif.stars=TRUE, cs.ind=1:2, tst.ind=3, P.values=TRUE, has.Pvalue=TRUE)
      cat("\n")
    }
  }
  cat("\n\n")

  invisible(x)
}

predBeta <- function(obj, p){
  fun <- obj$internal$fun
  q <- obj$internal$q
  Beta <- list()
  B <- matrix(0, nrow=length(p), ncol=5, dimnames=list(1:length(p), c("p", "beta", "se", "low", "up")))
  B[, 1] <- p
  bet <- fun(obj$x, p)

  gr <- matrix(0, length(p), length(obj$x))
  gr <- grad2(fun, x=obj$x, p=p, n.var=obj$internal$n.var, q=q, s=obj$s, ind=obj$ind, h=obj$internal$h1, method=obj$internal$meth)$grad
  if(length(p) == 1) gr <- t(gr)
  s2 <- c(0, cumsum(obj$s))
  lowBet <- upBet <- seBet <- matrix(0, nrow=length(p), ncol=q, dimnames=list(p, c("(Intercept)", paste0("X", 1:(q-1)))))

  for(j in 1:length(p)){
    for(i in 1:q){
      seBet[j, i] <- sqrt(t(gr[j, (s2[i]+1):s2[i+1]]) %*% obj$covar[(s2[i]+1):s2[i+1], (s2[i]+1):s2[i+1]] %*% as.matrix(gr[j, (s2[i]+1):s2[i+1]]))
      lowBet[j, i] <- bet[j, i]-1.96*seBet[j, i]
      upBet[j, i] <- bet[j, i]+1.96*seBet[j, i]
    }
  }

  for(i in 1:q){
    B[, 2] <- bet[, i]
    B[, 3] <- seBet[, i]
    B[, 4] <- lowBet[, i]
    B[, 5] <- upBet[, i]
    Beta[[i]] <- B
  }

  names(Beta) <- c("(Intercept)", paste0("X", 1:(q-1)))
  return(Beta)
}

#' @export
predict.niqr <- function(object, type = c("beta", "CDF", "QF", "sim"), newdata, p, ...){
  if(is.na(match(type <- type[1], c("beta", "CDF", "QF", "sim")))){stop("invalid 'type'")}

  if(type == "beta"){
    if(missing(p)){p <- seq.int(0.01,0.99,0.01)}
    if(any(p <= 0 | p >= 1)){stop("0 < p < 1 is required")}
    return(predBeta(object, p))
  }

  if(type == "CDF"){
    # q <- ncol(object$internal$X)-1
    nomi <- colnames(object$internal$X)[-1]
    if(!missing(newdata)){
      newdata <- as.matrix(newdata)
      y <- newdata[,"y"]
      X <- newdata[, nomi]
      X <- cbind(1, X)
    }else{
      X <- object$internal$X
      y <- object$internal$y
    }
    return(CDF=bisec(object$internal$fun, object$internal$fun2, object$x, X, y, nit=15, n=nrow(X)))
  }

  if(type == "QF"){
    # q <- ncol(object$internal$X)-1
    nomi <- colnames(object$internal$X)[-1]
    if(!missing(newdata)){
      newdata <- as.matrix(newdata)
      y <- newdata[, "y"]
      X <- newdata[, nomi]
      X <- cbind(1, X)
    }else{
      X <- object$internal$X
      y <- object$internal$y
    }
    pp <- bisec(object$internal$fun, object$internal$fun2, object$x, X, y, nit=15, n=nrow(X))
    beta <- object$internal$fun(object$x, pp)
    return(QF=object$internal$fun2(beta, X=X))
  }

  if(type == "sim"){
    n <- object$internal$n
    beta <- object$internal$fun(object$x, runif(n))
    return(ysim=object$internal$fun2(beta, X=object$internal$X))
  }
}

#' @export
plot.niqr <- function(x, conf.int=TRUE, which=NULL, ask=TRUE, ...){

  plot.niqr.int <- function(p,u,j,conf.int,L){
    beta <- u[[j]][, 2]
    low <- u[[j]][, 4]
    up <- u[[j]][, 5]
    if(is.null(L$ylim)){
      if(conf.int){y1 <- min(low); y2 <- max(up)}
      else{y1 <- min(beta); y2 <- max(beta)}
      L$ylim <- c(y1,y2)
    }
    plot(p, beta, xlab = L$xlab, ylab = L$ylab[j], main = L$labels[j],
         type = "l", lwd = L$lwd, xlim = L$xlim, ylim = L$ylim, col = L$col, cex.lab = L$cex.lab, cex.axis = L$cex.axis)
    if(conf.int){
      points(p, low, lty = 2, lwd = L$lwd, type = "l", col = L$col, cex.lab = L$cex.lab, cex.axis = L$cex.axis)
      points(p, up, lty = 2, lwd = L$lwd, type = "l", col = L$col, cex.lab = L$cex.lab, cex.axis = L$cex.axis)
    }
  }

  L <- list(...)
  if(is.null(L$xlim)){L$xlim = c(0.01,0.99)}
  if(is.null(L$lwd)){L$lwd <- 2}
  if(is.null(L$cex.lab)){L$cex.lab <- 1}
  if(is.null(L$cex.axis)){L$cex.axis <- 1}
  if(is.null(L$col)){L$col <- "black"}
  if(is.null(L$xlab)){L$xlab <- "p"}
  if(is.null(L$labels)) L$labels <- colnames(x$internal$X)
  q <- length(L$labels)
  if(is.null(L$ylab)){L$ylab <- rep("beta(p)", q)}
  L$labels <- c(L$labels, "qqplot")

  p <- seq.int(max(0.001,L$xlim[1]), min(0.999,L$xlim[2]), length.out = 100)
  u <- predict(x, p=p, type="beta")

  if(!is.null(which) | !ask){
    if(is.null(which)){which <- 1:q}
    nf <- n2mfrow(length(which))
    par(mfrow=nf)
    for(j in which){plot.niqr.int(p,u,j,conf.int,L)}
    par(mfrow=c(1,1))
  }
  else{
    pick <- 1
    while(pick > 0 && pick <= q + 1){
      pick <- menu(L$labels, title = "Make a plot selection (or 0 to exit):\n")
      if(pick > 0 && pick <= q){plot.niqr.int(p,u,pick,conf.int,L)}
      else if(pick == q + 1){
        KM <- NULL
        n <- x$internal$n
        KM$time <- 1:n/(n+1)
        KM$cdf <- sort(predict(x, type="CDF"))
        plot(KM$time, KM$cdf, pch = 20, cex = 0.5,
             xlim = c(0,1), ylim = c(0,1),
             ylab = "U(0,1) quantiles", xlab = "fitted CDF quantiles", cex.axis = L$cex.axis, cex.lab = L$cex.lab)
        abline(0,1)
      }
    }
  }
}

#' @export
test.fit.niqr <- function(obj, R = 100){

  test.unif.ct <- function(y){
    n <- length(y)
    o <- order(y)
    y <- y[o]
    w <- rep.int(1, n)
    W <- cumsum(w)
    hat.Fy <- W/W[n]
    Fy <- y
    DD <- Fy - hat.Fy
    ks <- max(abs(DD))

    Fy <- c(0, Fy, 1)
    hat.Fy <- c(0, hat.Fy, 1)
    y <- c(0, y, 1)
    U <- (hat.Fy - Fy)^2
    n <- n + 2
    h <- y[2:n] - y[1:(n-1)]
    b1 <- U[1:(n-1)]
    b2 <- U[2:n]
    A <- (b1 + b2)*h/2
    cvm <- sum(A)

    return(list(ks=ks, cvm=cvm))
  }

  n <- obj$internal$n
  q <- obj$internal$q
  y <- obj$internal$y
  X <- obj$internal$X
  p.star.y <- obj$p.star.y
  theta0 <- obj$x
  fun <- obj$internal$fun
  fun2 <- obj$internal$fun2

  test0 <- unlist(test.unif.ct(y = p.star.y))
  test <- matrix(NA, R, 2)

  pb <- txtProgressBar(min=0, max=R, style=3)
  for(b in 1:R){
    setTxtProgressBar(pb, b)

    id <- sample.int(n, size=n , replace=TRUE)
    Xb <- X[id,, drop = FALSE]
    beta <- fun(theta0, runif(n))
    yb <- fun2(beta, X=Xb)

    mod <- niqr(fun, fun2, obj$internal$x0, X=Xb, y=yb, control=list(tol=obj$internal$tol,
                alpha=obj$internal$alpha, beta=obj$internal$beta, maxit=obj$internal$maxit,
                maxit_start=obj$internal$maxit_start, low_p=obj$internal$pp[1],
                up_p=obj$internal$pp[obj$internal$npp], n_p=obj$internal$npp,
                h1=obj$internal$h1, meth=obj$internal$meth, display=FALSE))
    p.star.y.new <- mod$p.star.y

    test[b, ] <- unlist(test.unif.ct(y=p.star.y.new))
  }
  close(pb)

  out <- cbind(test0*c(1, n), c(mean(test[,1] >= test0[1], na.rm = TRUE), mean(test[,2] >= test0[2], na.rm = TRUE)))
  rownames(out) <- c("Kolmogorov-Smirnov", "Cramer-Von Mises")
  colnames(out) <- c("statistic", "p-value")

  out
}

grad2 <- function(fun, x, p, n.var, q, s, ind, h=1e-4, method=c("1", "2", "3")){
  # ff <- function(x, p=p, X=X, fun=fun) fun2(x, p, X, fun)
  method <- match.arg(method)

  grad <- switch(method,
         "1"={grad <- NULL
              f0 <- fun(x, p)
              for(i in 1:n.var) grad <- cbind(grad, (fun(x + c(rep(0, (i-1)), h, rep(0, n.var-i)), p) - f0)/h)
              grad},
         "2"={grad <- NULL
              for(i in 1:n.var){
                ff <- c(rep(0, (i-1)), h, rep(0, n.var-i))
                grad <- cbind(grad, (fun(x + ff, p) - fun(x - ff, p))/(2*h))}
              grad},
         "3"={grad <- NULL
              for(i in 1:n.var){
                ff <- c(rep(0, (i-1)), h, rep(0, n.var-i))
                grad <- cbind(grad, (-fun(x + 2*ff, p) + 8*fun(x + ff, p) - 8*fun(x - ff, p) + fun(x - 2*ff, p))/(12*h))}
              grad})

  if(is.null(s) | is.null(ind)){
    grad_2 <- grad
    colnames(grad_2) <- 1:dim(grad)[2]
    ind <- as.integer(names(unlist(apply(grad_2, 2, function(.x) which(sum(abs(.x)) > 1e-9)))))
    s <- table(colnames(grad[, ind]))
  }

  grad <- grad[, ind]

  return(list(grad=grad, s=s, ind=ind))
}

grad3 <- function(fun, fun2, x, p, X, npp, q, h=1e-4, method=c("1", "2", "3")){
  # ff <- function(x, p=p, X=X, fun=fun) fun2(x, p, X, fun)
  method <- match.arg(method)
  beta <- fun(x, p)

  grad <- switch(method,
                 "1"={grad <- NULL
                      f0 <- fun2(beta, X)
                      for(i in 1:q){
                        ff <- matrix(c(rep(0, (i-1)), h, rep(0, q-i)), ncol=q, nrow=npp, byrow=T)
                        grad <- cbind(grad, (fun2(beta + ff, X) - f0)/h)}
                      grad},
                 "2"={grad <- NULL
                      for(i in 1:q){
                        ff <- matrix(c(rep(0, (i-1)), h, rep(0, q-i)), ncol=q, nrow=npp, byrow=T)
                        grad <- cbind(grad, (fun2(beta + ff, X) - fun2(beta - ff, X))/(2*h))}
                      grad},
                 "3"={grad <- NULL
                      for(i in 1:q){
                        ff <- matrix(c(rep(0, (i-1)), h, rep(0, q-i)), ncol=q, nrow=npp, byrow=T)
                        grad <- cbind(grad, (-fun2(beta + 2*ff, X) + 8*fun2(beta + ff, X) - 8*fun2(beta - ff, X) + fun2(beta - 2*ff, X))/(12*h))}
                      grad})
  # grad <- grad[, which(abs(grad[1,]) >= 1e-9)]

  return(grad)
}

num.int <- function(x, dx, fx){
  n <- length(dx) + 1
  k <- ncol(fx)
  fL <- fx[1:(n-1),, drop = FALSE]
  fR <- fx[2:n,, drop = FALSE]

  out <- apply(rbind(0, 0.5*dx*(fL + fR)), 2, cumsum)

  OUT <- list()
  for(j in 1:k){
    fL <- out[2:n,j]; fR <- out[1:(n-1),j]
    if(all(fL <= fR) | all(fL >= fR)){
      method <- "hyman"
    }else{
      method <- "fmm"
    }
    OUT[[j]] <- splinefun(x, out[,j], method = method)
  }
  OUT
}

integr <- function(fun, a=.01, b=.99, n=100, p0, h, ...){
  f <- match.fun(fun)
  ff <- function(x) f(x, ...)

  if(missing(p0)){
    h <- (b-a)/n
    seqK <- seq(n-1)
    p0 <- c(a, a + seqK*h, b)
  }else{
    n <- length(p0)
    a <- p0[1]
    b <- p0[n]
    if(missing(h)) h <- (b - a)/n
  }
  int <- h * c(ff(a)/2, ff(p0[-c(1, n)]), ff(b)/2)
  sum(int[which(is.finite(int))])

  # if(!is.finite(ff(a))) a <- a + 1e-4
  # if(!is.finite(ff(b))) b <- b - 1e-4
  # h * (ff(a)/2 + sum(ff(p0)) + ff(b)/2)
}

bisec <- function(fun, fun2, x, X, y, nit=15, n=NULL){
  p0 <- rep(.5, n)
  for(i in 1:nit){
    beta <- fun(x, p0)
    p1 <- p0 + sign(y - fun2(beta, X)) * (.5^(i+1))
    p0 <- p1
  }
  return(p1)
}

Bfun <- function(fun, fun2, x, p.star.y, xx, dx, npp, q, X, n.var, p0, s=NULL, ind=NULL, hh, a=.01, b=.99, n=100,
                 fun_prime_theta=NULL, fun_prime_beta=NULL, ...){
    if(missing(xx)) xx <- seq(a, b, length=n)
    if(missing(dx)) dx <- xx[2:n] - xx[1:(n-1)]

    if(is.null(fun_prime_beta)){
      dhdb <- grad3(fun, fun2, x, xx, X, npp, q, ...)
    }else{
      beta <- fun(x, xx)
      dhdb <- fun_prime_beta(beta, X)
    }
    if(is.null(fun_prime_theta)){
      dbdt <- grad2(fun, x, xx, n.var, q, s, ind, ...)
      s <- dbdt$s
      ind <- dbdt$ind
      dbdt <- dbdt$grad
    }else{
      dbdt <- fun_prime_theta(x, xx)
      ind <-  NULL
      s <- dbdt[[2]]
      dbdt <- dbdt[[1]]
    }
    s2 <- c(0, cumsum(s))

    A <- num.int(x=xx, dx=dx, fx=dbdt)
    BB <- sapply(seq(length(A)), function(i) integr(function(p, i) A[[i]](p), a=a, b=b, n=n, p0=p0, h=hh, i=i))

    S1 <- NULL
    for(i in 1:q) for(j in (s2[i]+1):s2[i+1]){S1 <- cbind(S1, dhdb[, i]*(BB[j]-A[[j]](p.star.y)))}

    return(list(S1=S1, s=s, ind=ind, dhdb=dhdb, dbdt=dbdt))
} #fun = funLin, fun2=Q

bfun <- function(fun, fun2, x, p.star.y, xx, dx, npp, q, X, p0, hh, a=.01, b=.99, n=100){
    if(missing(xx)) xx <- seq(a, b, length=n)
    if(missing(dx)) dx <- xx[2:npp] - xx[1:(npp-1)]
    beta <- fun(x, xx)
    fxx <- fun2(beta, X)
    A <- num.int(x=xx, dx=dx, fx=fxx)

    BB <- Bh <- NULL
    for(i in 1:q){
      Bh <- cbind(Bh, A[[i]](p.star.y))
      BB <- c(BB, integr(function(p, i) A[[i]](p), a=a, b=b, n=n, p0=p0, h=hh, i=i))
    }

    S1 <- t(apply(Bh, 1, function(.x) BB - .x))

    return(S1=S1)
}

funGrad <- function(fun, fun2, x, y, X, i=FALSE, p.star.y, xx=NULL, dx=NULL, npp=NULL, p0=NULL, s=NULL, ind=NULL,
                    h=NULL, n=NULL, q=NULL, n.var=NULL, h1=1e-2, method=c("1","2","3"), fun_prime_theta=NULL,
                    fun_prime_beta=NULL){
  method <- match.arg(method)
  if(missing(p.star.y)) p.star.y <- bisec(fun, fun2, x, X, y, n=n)
  step <- Bfun(fun, fun2, p.star.y=p.star.y, x=x, xx=xx, dx=dx, npp=npp, q=q, X=X, n.var=n.var, p0=p0, s=s, ind=ind, hh=h, h=h1, method=method,
               fun_prime_theta=fun_prime_theta, a=xx[1], b=xx[npp], n=npp)
  s <- step$s
  ind <- step$ind
  S1 <- step$S1

  if(!i){
    grad <- colMeans(S1)
  }else{
    grad <- S1/n
  }

  # B1 <- sapply(seq(k), function(i) A[[i]](p.star.y, deriv=1))
  # prova <- NULL; for(j in 1:n.var) prova <- cbind(prova, rowSums(X * B1[, ind_i[j]:ind_f[j]]))
  # A1 <- 1/rowSums(prova)
  # Xw <- X*A1
  # b <- fun.p.theta(x, p.star.y)
  # J <- NULL
  # for (i1 in 1:n.var) {
  #   h.temp <- NULL
  #   for (i2 in 1:n.var) {
  #     h.temp <- cbind(h.temp, sum(crossprod(Xw, (X * (Bh[, ind_i[i2]:ind_f[i2]] * Bh[, ind_i[i1]:ind_f[i1]])))))
  #   }
  #   J <- rbind(J, h.temp)
  # }

  return(list(grad=grad, s=s, ind=ind))#, J=J))
}

Fun <- function(fun, fun2, x, y, X, q, n, p0, npp){
  f <- function(p, x, y, X, q, N){beta <- fun(x, p); qq <- y - fun2(beta, X); sum((p - (qq < 0))*qq)/N}
  f0 <- 0
  for(i in 1:npp) f0 <- f0 + f(p0[i], x, y, X, q, n)
  f0/npp
}

armijo <- function(x0, d, dtg, t0, f0, fun, fun2, alpha, beta, X, y, xx, dx, p0, h, n, npp, q){
  iter <- 100
  t <- t0
  xn <- x0 + t*d
  # p.star.y <- bisec(fun, fun2, xn, X=X, y=y, n=n, q=q)
  # fn <- Fun(fun, fun2, xn, X=X, y=y, xx=xx, dx=dx, npp=npp, p0=p0, h=h, n=n, q=q)
  fn <- Fun(fun, fun2, xn, y, X, q, n, p0, npp)
  if(!is.finite(fn)) fn <- f0
  i <- 0
  # g <- NULL
  while((fn > f0 + alpha*t*dtg) & (i <= iter)){
    t <- beta*t
    xn <- x0 + t*d
    # p.star.y <- bisec(fun, fun2, xn, X=X, y=y, n=n)
    # fn <- Fun(fun, fun2, xn, X=X, y=y, xx=xx, dx=dx, npp=npp, p0=p0, h=h, n=n, q=q)
    fn <- Fun(fun, fun2, xn, y, X, q, n, p0, npp)
    # g <- as.matrix(funGrad(fun, fun2, xn, y, X, i=FALSE, p.star.y, xx, dx, npp, p0, s, ind, h, n, q, n.var, h1, meth)$grad)
    # descent <- - H %*% g
    # DTG <- sum(descent*g)
    # const <- (abs(DTG) > beta * abs(dtg))
    if(!is.finite(fn)) fn <- f0
    i <- i + 1
  }
  # if(is.null(g)){
  #   p.star.y <- bisec(fun, fun2, xn, X=X, y=y, n=n)
  #   fn <- Fun(fun, fun2, xn, y, X, q, n, p0, npp)
  #   g <- as.matrix(funGrad(fun, fun2, xn, y, X, i=FALSE, p.star.y, xx, dx, npp, p0, s, ind, h, n, q, n.var, h1, meth)$grad)
  # }

  # t <- seq(1e-6, 1, l=100)
  # xn <- sapply(seq(100), function(i) x0 + t[i]*d)
  # fn <- apply(xn, 2, function(.x) Fun(fun, fun2, .x, y, X, q, n, p0, npp))
  # ind <- which.min(fn)

  obj <- list(t=t, xnew=xn, fn=fn)#, gn=g, p.star.y=p.star.y)
  return(obj)
}

start.theta.niqr <- function(fun, fun2, x0, X, y, control=list()){
  con <- list(tol=1e-6, alpha=.1, beta=.5, maxit=200, maxit_start=1, cluster=NULL, display=FALSE, epsilon=1e-12,
              a1=.001, h1=1e-2, meth="3", low_p=.0001, up_p=.9999, n_p=100, fun_prime_theta=NULL, fun_prime_beta=NULL)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms, collapse = ", "))

  if(missing(fun2) & is.null(con$fun_prime_beta)) con$fun_prime_beta <- funQ1
  if(missing(fun2)) fun2 <- Q
  h1 <- con$h1
  meth <- con$meth
  alpha <- con$alpha
  beta <- con$beta
  maxit_start <- con$maxit_start
  X <- as.matrix(X)
  q <- dim(X)[2]
  n <- dim(X)[1]
  if(ncol(X) == 1){
    colnames(X) = "(Intercept)"
  }else{
    if(is.null(colnames(X[, -1]))){
      colnames(X) <- c("(Intercept)", paste0("X", 1:(q-1)))
    }else{
      colnames(X)[1] <- "(Intercept)"
    }
  }
  # p10 <- seq.int(1/1024, 1023/1024, length.out = 1023) #c(0.1^(6:4), p10, 1 - 0.1^(4:6)) #seq(.01, .99, l=100) #(1:49/50)
  xx <- pp <- seq(con$low_p, con$up_p, l=con$n_p)
  npp <- con$n_p
  dx <- xx[2:npp] - xx[1:(npp-1)]
  h <- (pp[npp]-pp[1])/npp
  seqK <- seq(npp-1)
  p0 <- c(pp[1], pp[1]+seqK*h, pp[npp])

  n.var <- length(x0)
  x <- x0

  p.star.y <- bisec(fun, fun2, x, X=X, y=y, n=n)
  # f <- Fun(fun, fun2, x, X=X, y=y, p.star.y=p.star.y, xx=xx, dx=dx, npp=npp, p0=p0, h=h, n=n, q=q)
  f <- Fun(fun, fun2, x, y, X, q, n, p0, npp)
  gr <- funGrad(fun, fun2, x, y=y, X=X, p.star.y=p.star.y, xx=xx, dx=dx, npp=npp, p0=p0,
                s=NULL, ind=NULL, h=h, n=n, q=q, n.var=n.var, h1=h1, method=meth,
                fun_prime_theta=con$fun_prime_theta, fun_prime_beta=con$fun_prime_beta)
  s <- gr$s
  ind <- gr$ind
  g <- as.matrix(gr$grad)

  nmg <- max(abs(g))
  it <- 0
  a1 <- con$a1
  epsilon <- con$epsilon
  tol <- con$tol
  display <- con$display
  while(nmg > tol){ # while(nmg > tol*min(1, nmg0)){
    it <- it+1
    if(it > maxit_start){break}
    if(display) cat("Iter = ",formatC(it, digits=0, width=4, format="d"),
                    "  f = ",formatC(f, digits=8, width=11, format="f"),
                    "  sup(g) = ",formatC(nmg, digits=12, width=14, format="f"),"\n",sep="")
    if(it <= 10){
      H <- diag(1, nrow=n.var)
    }else{
      Dx <- x - xpr # sk
      Dg <- g - gpr # yk

      Dxg <- c(crossprod(Dx, Dg))
      DgHdg <- c(crossprod(Dg, H) %*% Dg)

      if((Dxg < epsilon) | (DgHdg < epsilon)){
        H <- diag(1, nrow=n.var)
      }else{
        HDg <- H %*% Dg
        dxdx <- tcrossprod(Dx, Dx)
        HdgHdg <- tcrossprod(HDg, HDg)

        u <- Dx / Dxg - HDg / DgHdg
        uu <- tcrossprod(u)

        H <- H + dxdx/Dxg - HdgHdg/DgHdg + DgHdg * uu
      }
    }
    descent <- -H %*% g   #pk
    dtg <- sum(descent*g) #grad(f)^T * pk
    # if(dtg > -min(a1, nmg)*nmg*sum(descent^2)){
    #   descent <- -g
    #   dtg <- sum(descent*g)
    # }
    obj <- armijo(x, descent, dtg, 1, f, fun, fun2, alpha, beta, X, y, xx, dx, p0, h, n, npp, q)
    xpr <- x
    x <- obj$xn

    gpr <- g
    p.star.y <- obj$p.star.y
    p.star.y <- bisec(fun, fun2, x, X=X, y=y, n=n)
    f <- obj$fn
    g <- as.matrix(funGrad(fun, fun2, x, y=y, X=X, p.star.y=p.star.y, xx=xx, dx=dx, npp=npp, p0=p0, s=s, ind=ind, h=h, n=n, q=q, n.var=n.var, h1=h1, method=meth,
                           fun_prime_theta=con$fun_prime_theta, fun_prime_beta=con$fun_prime_beta)$grad)
    nmg_temp <- max(abs(g))

    if(abs(nmg - nmg_temp) < epsilon) break
    nmg <- nmg_temp
  }

  objj <- list(x=x, fx=f, gx=g, Hx=H, nmg=nmg, p.star.y=p.star.y, s=s, ind=ind, iter=it, call=call)
  objj$internal <- list(q=q, n=n, xx=xx, pp=pp, npp=npp, dx=dx, h=h, seqK=seqK, p0=p0, n.var=n.var, y=y, X=X, h1=h1, meth=meth, epsilon=epsilon)

  return(objj)
}

#' @export
niqr <- function(fun, fun2, x0,  X, y, control=list()){
  call <- match.call()

  con <- list(tol=1e-6, alpha=.1, beta=.5, maxit=200, maxit_start=1, cluster=NULL, display=FALSE, epsilon=1e-12,
              a1=.001, h1=1e-2, meth="3", low_p=.0001, up_p=.9999, n_p=100, fun_prime_theta=NULL, fun_prime_beta=NULL)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms, collapse = ", "))

  if(missing(fun2) & is.null(con$fun_prime_beta)) con$fun_prime_beta <- funQ1
  if(missing(fun2)) fun2 <- Q
  ob <- start.theta.niqr(fun, fun2, x0, X=X, y=y, control=con)

  maxit <- con$maxit
  maxit_start <- con$maxit_start
  a1 <- con$a1
  epsilon <- con$epsilon
  tol <- con$tol
  npp <- con$n_p
  h1 <- con$h1
  meth <- con$meth
  display <- con$display
  alpha <- con$alpha
  beta <- con$beta
  cluster <- con$cluster

  q <- ob$internal$q
  n <- ob$internal$n
  xx <- pp <- ob$internal$pp
  dx <- ob$internal$dx
  h <- ob$internal$h
  seqK <- ob$internal$seqK
  p0 <- ob$internal$p0
  X <- ob$internal$X
  y <- ob$internal$y
  n.var <- ob$internal$n.var

  x <- ob$x
  p.star.y <- ob$p.star.y
  f <- ob$fx
  g <- ob$gx
  nmg <- ob$nmg
  s <- ob$s
  ind <- ob$ind

  it <- 0
  code <- 0
  while(nmg > tol){ # while(nmg > tol*min(1, nmg0)){
    it <- it+1
    if(it > maxit){code <- 1; break}
    if(display) cat("Iter = ",formatC(it, digits=0, width=4, format="d"),
                    "  f = ",formatC(f, digits=8, width=11, format="f"),
                    "  sup(g) = ",formatC(nmg, digits=12, width=14, format="f"),"\n",sep="")
    if(it == 1){
      H <- diag(diag(ob$Hx), nrow=n.var)
      descent <- - H %*% g
      dtg <- sum(descent*g)
    }else{
      Dx <- x - xpr
      Dg <- g - gpr

      Dxg <- c(crossprod(Dx, Dg))
      DgHdg <- c(crossprod(Dg, H) %*% Dg)

      if((Dxg < epsilon) | (DgHdg < epsilon)){
        descent <- -g
        dtg <- sum(descent*g)
      }else{
        HDg <- H %*% Dg
        dxdx <- tcrossprod(Dx, Dx)
        HdgHdg <- tcrossprod(HDg, HDg)

        u <- Dx / Dxg - HDg / DgHdg
        uu <- tcrossprod(u)

        H <- H + dxdx/Dxg - HdgHdg/DgHdg + DgHdg * uu

        # H <- H + tcrossprod(Dg)/c(crossprod(Dg, Dx)) - tcrossprod(H %*% Dx)/c(crossprod(Dx, H) %*% Dx) #BFGS
        # H <- H + tcrossprod(Dx)/c(crossprod(Dx, Dg)) + tcrossprod(H %*% Dg)/c(crossprod(Dg, H) %*% Dg) #DFP2
        # H <- H + (1 + DgHdg/Dxg) * tcrossprod(Dx)/Dxg - (Hdgdx + t(Hdgdx))/Dxg

        descent <- -H %*% g
        dtg <- sum(descent*g)
      }
    }
    # if(dtg > -min(a1, nmg)*nmg*sum(descent^2)){
    #   descent <- -g
    #   dtg <- sum(descent*g)
    # }
    obj <- armijo(x, descent, dtg, 1, f, fun, fun2, alpha, beta, X, y, xx, dx, p0, h, n, npp, q)
    xpr <- x
    x <- obj$xn

    gpr <- g
    p.star.y <- obj$p.star.y
    p.star.y <- bisec(fun, fun2, x, X=X, y=y, n=n)
    f <- obj$fn
    g <- as.matrix(funGrad(fun, fun2, x, y=y, X=X, p.star.y=p.star.y, xx=xx, dx=dx, npp=npp, p0=p0, s=s, ind=ind, h=h, n=n, q=q, n.var=n.var, h1=h1, method=meth,
                           fun_prime_theta=con$fun_prime_theta, fun_prime_beta=con$fun_prime_beta)$grad)
    nmg_temp <- max(abs(g))

    # g <- as.matrix(funGrad(fun, fun2, x, y=y, X=X, p.star.y=p.star.y, xx=xx, dx=dx, npp=npp, p0=p0, s=s, ind=ind, h=h, n=n, q=q, n.var=n.var, h1=h1, method=meth)$grad)

    if(abs(nmg - nmg_temp) < epsilon){code <- 2; break}
    nmg <- nmg_temp
  }

  if(it <= 10) H <- ob$Hx

  # func <- function(x, y, X, q, n, p0, npp){
  #   f <- function(p, x, y, X, q, N){beta <- fun(x, p); qq <- y - fun2(beta, X); sum((p - (qq < 0))*qq)/N}
  #   f0 <- 0
  #   for(i in 1:npp) f0 <- f0 + f(p0[i], x, y, X, q, n)
  #   f0/npp
  # }
  # H <- solve(rootSolve::hessian(func, x=x, centered=TRUE, pert=1e-2, y=y, X=X, q=q, n=n, p0=p0, npp=npp))

  if(is.null(cluster)){
    matGrad <- funGrad(fun, fun2, x, y=y, X=X, p.star.y=p.star.y, i=TRUE, xx=xx, dx=dx, npp=npp, p0=p0, s=s, ind=ind, h=h, n=n, q=q, n.var=n.var, h1=h1, method=meth,
                       fun_prime_theta=con$fun_prime_theta, fun_prime_beta=con$fun_prime_beta)$grad
    # s <- matGrad$s
    # matGrad <- matGrad$grad
    omega <- (t(matGrad) %*% matGrad)
    covar <- H %*% omega %*% H
    se <- sqrt(diag(covar))
  }else{
    M <- length(unique(cluster))
    dfc <- (M/(M-1))*((n-1)/(n-n.var))
    matGrad <- funGrad(fun, fun2, x, y=y, X=X, p.star.y=p.star.y, i=TRUE, xx=xx, dx=dx, npp=npp, p0=p0, s=s, ind=ind, h=h, n=n, q=q, n.var=n.var, h1=h1, method=meth,
                       fun_prime_theta=con$fun_prime_theta, fun_prime_beta=con$fun_prime_beta)$grad
    # s <- matGrad$s
    # matGrad <- matGrad$grad
    matGrad <- apply(matGrad, 2, function(.x) tapply(.x, cluster, sum))
    omega <- (t(matGrad) %*% matGrad)
    covar <- dfc*(H %*% omega %*% H)
    se <- sqrt(diag(covar))
  }

  objj <- list(x=c(x), se=se, fx=f, gx=g, Hx=H, omega=omega, covar=covar, p.star.y=p.star.y, s=s, ind=ind, code=code, iter=it, call=call)
  objj$internal <- list(q=q, n=n, x0=x0, alpha=alpha, beta=beta, tol=tol, display=display, xx=xx, pp=pp, npp=npp, dx=dx, h=h, seqK=seqK,
                        p0=p0, n.var=n.var, y=y, X=X, h1=h1, meth=meth, maxit=maxit, maxit_start=maxit_start, fun=fun, fun2=fun2,
                        fun_prime_theta=con$fun_prime_theta, fun_prime_beta=con$fun_prime_beta)
  class(objj) <- "niqr"

  return(objj)
}
