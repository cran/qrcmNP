#' @importFrom stats integrate splinefun model.response model.weights model.matrix terms model.frame delete.response coef
#' @importFrom stats approxfun sd prcomp approx predict .getXlevels runif printCoefmat pnorm qexp pbeta qbeta pchisq
#' @importFrom survival Surv survfit coxph
#' @importFrom grDevices adjustcolor gray.colors
#' @importFrom graphics plot abline axis matplot text points
#' @importFrom utils getFromNamespace menu tail
#' @import qrcm

pmax0 <- function(x){(x + abs(x))/2}

#' @export
piqr <- function(formula, formula.p = ~ slp(p, 1), weights, data, s, nl=70,
                 display=TRUE, tol=1e-6, maxit=100){

  p.bisec.internal <- getFromNamespace("p.bisec.internal", "qrcm")
  start.theta <- getFromNamespace("start.theta", "qrcm")

  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "weights", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  A <- check.in.2(mf, formula.p, s)
  type <- A$type
  X <- A$U$X
  y <- A$U$y
  V <- A$V
  n <- nrow(X)
  p <- ncol(X)
  bfun <- A$internal.bfun
  attributes(A$bfun) <- c(attributes(A$bfun), A$stats.B)
  s <- A$s
  k <- ncol(s)
  intercept <- A$stats.X$intercept
  if(intercept) s_int <- matrix(as.integer(s[1, ]), nrow=1)
  if(missing(data)) data <- data.frame(X=X, y=y)

  covar <- seqS <- list()
  if(type == "iqr") minimum <- vector("double", length=nl) else minimum <- rep(NA, nl)
  beta <- matrix(0, nrow=p*k, ncol=nl)
  dl <- matrix(0, p*k, nl)
  flag <- TRUE

  theta0 <- start.theta(V$y, V$z, V$d, V$X, V$weights, bfun,
                        df = max(5, min(15, round(n/30/(q + 1)))), yy=NULL, zz=NULL, s=s)
  Theta0 <- matrix(0, ncol=k, nrow=p)

  if(!intercept){
    if(type == "iqr"){ee <- piqr.ee}
    else if(type == "ciqr"){ee <- pciqr.ee}
    else{ee <- pctiqr.ee}

    oIQR <- NULL
    if(type == "iqr") minimum[1] <- NaN
    seqS[[1]] <- theta <- matrix(0, nrow=p, ncol=k)
    covar[[1]] <- matrix(0, nrow=p*k, ncol=p*k)
    p.star.y <- p.bisec.internal(theta, V$y, V$X, bfun$bp)
    p.star.z <-  NULL
    grad <- ee(theta0, V$y, V$z, V$d, V$X, V$Xw, bfun, p.star.y, p.star.z, J=F, segno=1, lambda=0)$g
    grad[!s] <- 0
  }else{
    seqS[[1]] <- A$s <- matrix(c(s_int, rep(0, (p-1)*k)), ncol=k, byrow=T)
    Theta0[A$s==1] <- theta0[A$s==1]
    oIQR <- piqr.internal(mf=mf, cl=cl, formula.p=formula.p, tol=tol, maxit=maxit,
                          segno=1, lambda=0, check=FALSE, A=A, s=seqS[[1]], st.theta=FALSE, theta0=Theta0)
    if(type == "iqr") minimum[1] <- oIQR$obj.function
    beta[, 1] <- c(oIQR$coefficients)
    covar[[1]] <- oIQR$covar
    grad <- attributes(oIQR$mf)$grad
    grad[!s] <- 0
  }

  np <- 1

  dl[, np] <- grad
  segno <- ifelse(dl[, np] >= 0, 1, -1)

  l <- max(abs(dl[, np]))
  lmin <- 1e-4 * l
  lratio <- exp((log(l[1]) - log(lmin)) / (nl - 1))
  for(i in 2:nl) l <- append(l, l[i-1] / lratio)

  al <- l[np]
  AA <- ifelse(abs(dl[, 1]) >= al, TRUE, FALSE)
  if(intercept) AA[seq(1, (p*k), by=p)] <- as.logical(s_int)

  if(display) cat("\n Iter=",np," lambda=",round(al,8),sep="")
  for(kk in 1:10000){
    if(flag){
      # print(max(abs(dl[, np])) - al)
      # if((max(abs(dl[, np])) - al) >= 1e-1) break

      np <- np + 1
      if(np > nl){np <- np-1;break}
      al <- l[np]
      flag <- TRUE
      if(display) cat("\n Iter=",np," lambda=",round(al,8),sep="")
    }else{
      flag <- TRUE
    }

    seqS[[np]] <- A$s <- matrix(as.integer(AA), ncol=k)
    Theta0[A$s==1] <- theta0[A$s==1]
    tempSegno <- A$s*segno
    if(intercept) tempSegno[1, ] <- rep(0, k)

    oIQR <- try(piqr.internal(mf=mf, cl=cl, formula.p=formula.p, tol=tol, maxit=maxit,
                              segno=tempSegno, lambda=al, check=FALSE, A=A, s=seqS[[np]], st.theta=FALSE, theta0=Theta0), silent = TRUE)
    if(class(oIQR) == "try-error"){np <- np-1; break}

    if(type == "iqr") minimum[np] <- oIQR$obj.function
    beta[, np] <- c(oIQR$coefficients)
    covar[[np]] <- oIQR$covar
    grad <- attributes(oIQR$mf)$grad
    grad[!s] <- 0

    dl[, np] <- grad

    if(all(A$s == s) & (max(abs(dl[, np])) < (min(l)*10))) break

    for(m in 1:(p*k)){
      if(!AA[m]){
        if(abs(dl[m, np]) >= al){
          flag <- FALSE
          AA[m] <- TRUE
          segno[m] <- ifelse(dl[m, np] >= 0, 1, -1)
        }
      }
    }
  }

  attr(mf, "assign") <- A$stats.X$assign
  attr(mf, "stats") <- list(B = A$stats.B, X = A$stats.X, y = A$stats.y)
  attr(mf, "all.vars") <- A$V
  attr(mf, "bfun") <- A$bfun
  attr(mf, "type") <- type

  # df <- apply(beta[, 1:np], 2, function(.x) length(which(.x != 0)))
  df_max <- apply(s, 1, sum)
  df <- vector(length=np)
  se <- coefficients <- list()
  for(i in 1:np){
    dimnames(seqS[[i]]) <- list(A$stats.X$coef.names, A$stats.B$coef.names)
    coefficients[[i]] <- matrix(beta[, i], ncol=k, dimnames=list(A$stats.X$coef.names, A$stats.B$coef.names))
    # df[i] <- length(which(beta[, i] != 0))
    df[i] <- sum(apply(coefficients[[i]], 1, function(.x) sum(I(.x != 0)))/df_max)
  }

  fit <- list(call=cl, lambda=l[1:np], coefficients=coefficients, covar=covar,
              obj.function=minimum[1:np], dl=dl[, 1:np], df=round(df, 3), seqS=seqS)
  fit$internal <- list(mf=mf, formula=formula, data=data, formula.p=formula.p, nlambda=np, X=X, y=y, k=k, intercept=intercept, type=type)

  class(fit) <- "piqr"
  return(fit)
}

#' @export
gof.piqr <- function(object, method=c("AIC","BIC","GIC","GCV","NIC"), plot=TRUE, ...){
  cl <- match.call()

  method <- match.arg(method)

  lambda <- object$lambda
  nl <- object$internal$nlambda
  X <- object$internal$X
  y <- object$internal$y

  intercept <- object$internal$intercept

  n <- dim(X)[1]
  p <- dim(X)[2]
  k <- object$internal$k
  df <- object$df
  beta <- object$coefficients
  obj.function <- object$obj.function

  if(method == "NIC"){
    if(intercept){
      dfDen <- df - length(which(beta[[nl]][1, ] != 0))
      Den <- sapply(1:nl, function(.i) sum(abs(beta[[.i]][-1, ])))
    }else{
      dfDen <- df
      Den <- sapply(1:nl, function(.i) sum(abs(beta[[.i]])))
    }
  }

  curve <- switch(method,
                  "AIC"=log(obj.function) + df/n,
                  "BIC"=log(obj.function) + log(n)*df/(2*n),
                  "GIC"=log(obj.function) + log(log(n))*log(p)*df/(2*n),
                  "GCV"=log(obj.function)/((n - df)^2),
                  "NIC"=sqrt(obj.function)/(Den/dfDen))

  posMinLambda <- which.min(curve)
  minLambda <- lambda[posMinLambda]
  dfMinLambda <- df[posMinLambda]
  betaMin <- beta[[posMinLambda]]

  if(plot){
    color <- gray.colors(2)
    plot(log(lambda), curve, type="l", col=color[1], xlab=expression(log(lambda)),
         ylab=bquote(.(method)), ...)
    abline(v=log(minLambda), col=color[2], lty=2, ...)
    axis(3, at=log(minLambda), labels=paste("df=", dfMinLambda, sep=""), ...)
  }

  ogg <- list(call=cl, minLambda=minLambda, dfMinLambda=dfMinLambda,
              posMinLambda=posMinLambda, betaMin=betaMin)
  # class(ogg) <- "gof.piqr"
  return(ogg)
}

piqr.internal <- function(mf, cl, formula.p, tol=1e-6, maxit=100, s,
                          segno=1, lambda=0, check=FALSE, A, st.theta=FALSE, theta0){

  start.theta <- getFromNamespace("start.theta", "qrcm")
  iobjfun <- getFromNamespace("iobjfun", "qrcm")
  km <- getFromNamespace("km", "qrcm")
  p.bisec.internal <- getFromNamespace("p.bisec.internal", "qrcm")
  p.bisec <- getFromNamespace("p.bisec", "qrcm")
  apply_bfun <- getFromNamespace("apply_bfun", "qrcm")
  trans <- getFromNamespace("trans", "qrcm")

  n <- nrow(mf)
  if(check) A <- check.in.2(mf, formula.p, s)
  V <- A$V; U <- A$U; s <- A$s; type <- A$type
  mf <- A$mf; n <- nrow(mf)
  S <- list(B = A$stats.B, X = A$stats.X, y = A$stats.y)
  attributes(A$bfun) <- c(attributes(A$bfun), S$B)
  bfun <- A$internal.bfun

  if(missing(maxit)){maxit <- 10 + 10*sum(s)}
  else{maxit <- max(10, maxit)}

  if(st.theta){
    q <- length(S$X$vars)
    if(type != "iqr" | q > 0){
      Ty <- trans(V$z,V$y,V$d,V$weights, type)
      yy <- Ty$f(V$y)
      zz <- (if(type == "ctiqr") Ty$f(V$z) else V$z)
    }
    else{yy <- zz <- NULL}

    theta0 <- start.theta(V$y, V$z, V$d, V$X, V$weights, bfun,
                          df = max(5, min(15, round(n/30/(q + 1)))), yy, zz, s = s)
  }else{
    theta0 <- theta0
  }

  Fit <- NULL
  fit.ok <- FALSE
  safeit <- 5
  try.count <- 0
  eeTol <- 0.5
  eps0 <- 0.1

  while(!fit.ok){

    try.count <- try.count + 1

    fit <- piqr.newton(theta0, V$y, V$z, V$d, V$X, V$Xw,
                      bfun, s = s, type = type, tol = tol, maxit = maxit, safeit = safeit, eps0 = eps0, segno = segno, lambda = lambda)

    if(fit.ok <- (fit$rank == ncol(fit$jacobian) & max(abs(fit$ee)) < eeTol)){
      # covar <- try(cov.theta(fit$coefficients, V$y, V$z, V$d, V$X, V$Xw,
      #                        V$weights, bfun, fit$p.star.y, fit$p.star.z, type, s = s, segno = segno, lambda = lambda), silent = TRUE)
      # fit.ok <- (class(covar) != "try-error")
      fit.ok <- TRUE
    }

    # covar.ok <- (if(fit.ok){(qr(covar$Q)$rank == fit$rank)} else FALSE)
    if(fit.ok){break}
    # else if(fit.ok){Fit <- fit; Cov <- covar}

    # if(try.count > 10 && !is.null(Fit)){break}
    # if(try.count == 20){break}
    eeTol <- eeTol + 0.5
    safeit <- safeit + 2
    eps0 <- eps0/2
  }

  if(!fit.ok && is.null(Fit)){stop("unable to fit the model: this can be due to severe misspecification")}
  # if(!covar.ok){warning("the estimated covariance matrix is deemed to be singular")}
  if(!fit.ok){fit <- Fit}#; covar <- Cov}
  if(!fit$converged){warning("the algorithm did not converge")}

  # minimized loss function

  if(type == "iqr"){
    v <- (S$y$M - S$y$m)/10
    fit$obj.function <- iobjfun(fit$coef, V$y,V$X,V$weights, bfun, fit$p.star.y)*v
  }

  # fitted CDFs

  CDFs <- data.frame(CDF.y = fit$py,
                     CDF.z = (if(type == "ctiqr") fit$pz else NA))
  attr(CDFs, "km") <- km(CDFs$CDF.z, CDFs$CDF.y, V$d, V$weights, type)

  # output

  covar <- cov.theta.2(fit$coefficients, V$y, V$z, V$d, V$X, V$Xw,
            V$weights, bfun, fit$p.star.y, fit$p.star.z, type, s = s, segno = segno, lambda = lambda, J=fit$jacobian)

  attr(mf, "assign") <- S$X$assign
  attr(mf, "stats") <- S
  attr(mf, "CDFs") <- CDFs
  attr(mf, "all.vars") <- V
  attr(mf, "all.vars.unscaled") <- U
  attr(mf, "Q0") <- covar$Q0
  attr(mf, "internal.bfun") <- bfun
  attr(mf, "bfun") <- A$bfun
  attr(mf, "theta") <- fit$coefficients
  attr(mf, "type") <- type
  attr(mf, "grad") <- fit$grad
  attr(mf, "jacobian") <- fit$jacobian
  attr(mf, "type") <- type

  out <- check.out.2(fit$coefficients, S, covar = covar$Q)
  fit <- list(coefficients = out$theta, call = cl,
              converged = fit$converged, n.it = fit$n.it,
              obj.function = fit$obj.function,
              covar = out$covar, mf = mf, s = s)
  jnames <- c(sapply(attr(A$bfun, "coef.names"),
                     function(x,y){paste(x,y, sep = ":")}, y = S$X$coef.names))
  dimnames(fit$covar) <- list(jnames, jnames)
  dimnames(fit$coefficients) <- dimnames(fit$s) <- list(S$X$coef.names, S$B$coef.names)


  # CDF and PDF, precision ~ 1e-6

  fit$CDF <- p.bisec(fit$coef,U$y,U$X,A$bfun)
  b1 <- apply_bfun(A$bfun, fit$CDF, "b1fun")
  fit$PDF <- 1/c(rowSums((U$X%*%fit$coef)*b1))
  fit$PDF[attr(fit$CDF, "out")] <- 0
  attributes(fit$CDF) <- attributes(fit$PDF) <- list(names = rownames(mf))
  if(any(fit$PDF < 0)){warning("quantile crossing detected (PDF < 0 at some y)")}

  # finish

  class(fit) <- "iqr"
  fit
}

check.in.2 <- function(mf, formula.p, s){

  is.slp <- getFromNamespace("is.slp", "qrcm")
  make.bfun <- getFromNamespace("make.bfun", "qrcm")
  slp.basis <- getFromNamespace("slp.basis", "qrcm")
  apply_bfun <- getFromNamespace("apply_bfun", "qrcm")
  num.fun <- getFromNamespace("num.fun", "qrcm")

  if(!missing(s) && all(s == 0)){stop("'s' cannot be all zero")}
  explore.s <- function(s, dim){
    if(dim == 2){s <- t(s)}
    out <- 1
    if((r <- nrow(s)) > 1){
      for(j in 2:r){
        done <- FALSE; rj <- s[j,]
        for(h in 1:(j - 1)){
          if(all(rj == s[h,])){out[j] <- out[h]; done <- TRUE}
        }
        if(!done){out[j] <- max(out) + 1}
      }
    }
    out
  }

  # weights

  if(any((weights <- model.weights(mf)) < 0)){stop("negative 'weights'")}
  if(is.null(weights)){weights <- rep.int(1, nrow(mf)); alarm <- FALSE}
  else{
    alarm <- (weights == 0)
    sel <- which(!alarm)
    mf <- mf[sel,]
    weights <- weights[sel]
    weights <- weights/mean(weights)
  }
  if(any(alarm)){warning("observations with null weight will be dropped", call. = FALSE)}
  if((n <- nrow(mf)) == 0){stop("zero non-NA cases", call. = FALSE)}

    # y,z,d, weights

  zyd <- model.response(mf)
  type <- attributes(zyd)$type
  zyd <- cbind(zyd)

  if(is.null(type)){y <- zyd[,1]; z <- rep.int(-Inf,n); d <- rep.int(1,n); type <- fittype <- "iqr"}
  else if(type == "right"){
    y <- zyd[,1]; z <- rep.int(-Inf,n); d <- zyd[,2]
    type <- (if(any(d == 0)) "ciqr" else "iqr")
    fittype <- "ciqr"
  }
  else if(type == "counting"){
    z <- zyd[,1]; y <- zyd[,2]; d <- zyd[,3]; type <- fittype <- "ctiqr"
    if(all(z < min(y))){type <- (if(any(d == 0)) "ciqr" else "iqr")}
  }
  attr(type, "fittype") <- fittype
  if(!(any(d == 1))){stop("all data are censored")}

  # x and b(p)

  X <- model.matrix(attr(mf, "terms"), mf); q <- ncol(X)
  termlabelsX <- attr(attr(mf, "terms"), "term.labels")
  assignX <- attr(X, "assign")
  coefnamesX <- colnames(X)

  # p1 is used to evaluate the splinefuns. A non-evenly spaced grid, with more values on the tails.
  # p2 is for external use (p.bisec). A grid with p reachable by bisection on the p scale.
  # p3 is for internal use (p.bisec.internal). A grid with p reachable by bisection on the index scale.
  p1 <- pbeta(seq.int(qbeta(1e-6,2,2), qbeta(1 - 1e-6,2,2), length.out = 1000),2,2)
  p2 <- (1:1023)/1024
  p3 <- pbeta(seq.int(qbeta(1/(1000*n),2.5,2.5), qbeta(1 - 1/(1000*n),2.5,2.5), length.out = 1023),2.5,2.5)

  if((use.slp <- is.slp(formula.p))){
    k <- attr(use.slp, "k")
    intercept <- attr(use.slp, "intercept") 	# slp(0) = 0?
    intB <- attr(use.slp, "intB")			        # b(p) includes 1?
    assignB <- (1 - intB):k
    termlabelsB <- paste("slp", 1:k, sep = "")
    coefnamesB <- (if(intB) c("(Intercept)", termlabelsB) else termlabelsB)
    k <- k + intB
  }
  else{
    B <- model.matrix(formula.p, data = data.frame(p = c(p1,p2,p3)))
    B1 <- B[1:1000,, drop = FALSE]
    B2 <- B[1001:2023,, drop = FALSE]
    B3 <- B[2024:3046,, drop = FALSE]

    k <- ncol(B)
    assignB <- attr(B, "assign")
    termlabelsB <- attr(terms(formula.p), "term.labels")
    coefnamesB <- colnames(B)
  }
  if(missing(s)){s <- matrix(1,q,k)}
  else{
    if(any(dim(s) != c(q,k))){stop("wrong size of 's'")}
    if(any(s != 0 & s != 1)){stop("'s' can only contain 0 and 1")}
  }

  # x singularities (set s = 0 where singularities occur)
  # x is dropped as in a linear model, irrespective of s.

  vx <- qr(X); selx <- vx$pivot[1:vx$rank]
  if(vx$rank < q){s[-selx,] <- 0}

  # b(p) singularities. Dropped row by row, based on s

  if(!use.slp && qr(B3)$rank < k){
    u <- explore.s(s,1)
    for(j in unique(u)){
      sel <- which(s[which(u == j)[1],] == 1)
      if(length(sel) > 1){
        vbj <- qr(B3[,sel, drop = FALSE])
        if((rj <- vbj$rank) < length(sel)){
          s[u == j, sel[-vbj$pivot[1:rj]]] <- 0
        }
      }
    }
  }

  # location-scale statistics for x, b(p), and y

  ry <- range(y); my <- ry[1]; My <- ry[2]

  sX <- apply(X,2,sd); mX <- colMeans(X)
  intX <- (length((constX <- which(sX == 0 & mX != 0))) > 0)
  varsX <- which(sX > 0); zeroX <- which(sX == 0 & mX == 0)
  sX[constX] <- X[1,constX]; mX[constX] <- 0; sX[zeroX] <- 1
  if(length(constX) > 1){zeroX <- c(zeroX, constX[-1]); constX <- constX[1]}

  if(!use.slp){
    sB <- apply(B3,2,sd); mB <- colMeans(B3)
    intB <- (length((constB <- which(sB == 0 & mB != 0))) > 0); varsB <- which(sB > 0)
    if(length(varsB) == 0){stop("the quantile function must depend on p")}
    if(length(constB) > 1){stop("remove multiple constant functions from 'formula.p'")}
    if(any(sB == 0 & mB == 0)){stop("remove zero functions from 'formula.p'")}
    sB[constB] <- B3[1,constB]; mB[constB] <- 0
  }
  else{
    sB <- rep.int(1, k); mB <- rep.int(0, k)
    if(intB){constB <- 1; varsB <- 2:k}
    else{constB <- integer(0); varsB <- 1:k}
  }

  if(all(s[, varsB] == 0)){stop("the quantile function must depend on p (wrong specification of 's')")}
  if(!(theta00 <- ((intX & intB) && s[constX, constB] == 1)))
  {my <- 0; My <- sd(y)*5; mX <- rep.int(0,q)}
  else{for(j in varsX){if(any(s[j,] > s[constX,])){mX[j] <- 0}}}
  if(!intB | (intB && any(s[,constB] == 0))){mB <- rep.int(0,k)}

  # Create bfun (only used by post-estimation functions)

  if(!use.slp){
    bfun <- list()
    if(intB){bfun[[constB]] <- function(p, deriv = 0){rep.int(1 - deriv, length(p))}}
    for(j in varsB){bfun[[j]] <- make.bfun(p1,B1[,j])}
    names(bfun) <- coefnamesB
    attr(bfun, "k") <- k
  }
  else{
    bfun <- slp.basis(k - intB, intercept)
    if(!intB){bfun$a[1,1] <- bfun$A[1,1] <- bfun$AA[1,1] <- 0}
    attr(bfun, "intB") <- intB
    B2 <- apply_bfun(bfun,p2, "bfun")
  }

  attr(bfun, "bp") <- B2
  attr(bfun, "p") <- p2

  # first scaling of x, b(p), y

  U <- list(X = X, y = y, z = z)
  X <- scale(X, center = mX, scale = sX)
  y <- (y - my)/(My - my)*10
  z <- z0 <- (z - my)/(My - my)*10
  z[z < min(y)] <- -Inf
  if(!use.slp){B3 <- scale(B3, center = mB, scale = sB)}

  # principal component rotations that I can apply to x and b(p); second scaling

  rotX <- (diag(1,q))
  MX <- rep.int(0,q); SX <- rep.int(1,q)

  rotB <- (diag(1,k))
  MB <- rep.int(0,k); SB <- rep.int(1,k)

  # Create bfun

  p <- p3

  if(!use.slp){
    bp <- B3
    dp <- p[-1] - p[-1023]
    b1p <- num.fun(dp,bp, "der")
    Bp <- num.fun(dp,bp, "int")
    BBp <- num.fun(dp,Bp, "int")
    BB1 <- BBp[1023,]
  }
  else{
    k <- attr(bfun, "k")
    pp <- matrix(NA, 1023, k + 1)
    pp[,1] <- 1; pp[,2] <- p
    if(k > 1){for(j in 2:k){pp[,j + 1] <- pp[,j]*p}}
    bp <- tcrossprod(pp, t(bfun$a))
    b1p <- cbind(0, tcrossprod(pp[,1:k, drop = FALSE], t(bfun$a1[-1,-1, drop = FALSE])))
    pp <- cbind(pp, pp[,k + 1]*p, pp[,k + 1]*p^2)
    Bp <- tcrossprod(pp[,2:(k + 2)], t(bfun$A))
    BBp <- tcrossprod(pp[,3:(k + 3)], t(bfun$AA))
    BB1 <- colSums(bfun$AA)

    if(!intB){
      bp <- bp[,-1, drop = FALSE]
      b1p <- b1p[,-1, drop = FALSE]
      Bp <- Bp[,-1, drop = FALSE]
      BBp <- BBp[,-1, drop = FALSE]
      BB1 <- BB1[-1]
    }
  }
  BB1 <- matrix(rep(BB1, each = n), n)
  bpij <- NULL; for(i in 1:ncol(bp)){bpij <- cbind(bpij, bp*bp[,i])}

  internal.bfun <- list(p = p, bp = bp, b1p = b1p, Bp = Bp, BBp = BBp, BB1 = BB1, bpij = bpij)
  attr(internal.bfun, "pfun") <- approxfun(c(p[1], 0.5*(p[-1023] + p[-1])),p, method = "constant", rule = 2)

  # output. U = the original variables. V = the scaled/rotated variables.
  # stats.B, stats.X, stats.y = lists with the values use to scale/rotate

  stats.B <- list(m = mB, s = sB, M = MB, S = SB, rot = rotB, const = constB, vars = varsB,
                  intercept = intB, term.labels = termlabelsB, assign = assignB, coef.names = coefnamesB)
  stats.X <- list(m = mX, s = sX, M = MX, S = SX, rot = rotX, const = constX, vars = varsX,
                  intercept = intX, term.labels = termlabelsX, assign = assignX, coef.names = coefnamesX)
  stats.y <- list(m = my, M = My)

  V <- list(X = X, Xw = X*weights, y = y, z = z, d = d, weights = weights)
  if(type == "ctiqr"){V$z0 <- z0}
  list(mf = mf, U = U, V = V, stats.B = stats.B, stats.X = stats.X, stats.y = stats.y,
       internal.bfun = internal.bfun, bfun = bfun, s = s, type = type)
}

check.out.2 <- function(theta, S, covar){

  blockdiag <- function(A, d, type = 1){
    h <- nrow(A); g <- d/h
    if(type == 1){
      out <- diag(1,d)
      for(j in 1:g){ind <- (j*h - h  + 1):(j*h); out[ind,ind] <- A}
    }
    else{
      out <- matrix(0,d,d)
      for(i1 in 1:h){
        for(i2 in 1:h){
          ind1 <- (i1*g - g  + 1):(i1*g)
          ind2 <- (i2*g - g  + 1):(i2*g)
          out[ind1, ind2] <- diag(A[i1,i2],g)
        }
      }
      out <- t(out)
    }
    out
  }

  mydiag <- function(x){
    if(length(x) > 1){return(diag(x))}
    else{matrix(x,1,1)}
  }

  th <- cbind(c(theta))
  q <- nrow(theta)
  k <- ncol(theta)
  g <- q*k
  aX <- S$X; ay <- S$y; aB <- S$B
  cX <- aX$const; cB <- aB$const

  ##########################

  A <- blockdiag(mydiag(1/aX$S), g)
  th <- A%*%th
  covar <- A%*%covar%*%t(A)

  if(aX$intercept){
    A <- diag(1,q); A[cX,] <- -aX$M; A[cX, cX] <- 1
    A <- blockdiag(A,g)
    th <- A%*%th
    covar <- A%*%covar%*%t(A)
  }

  ##########################

  A <- blockdiag(mydiag(1/aB$S),g,2)
  th <- A%*%th
  covar <- A%*%covar%*%t(A)

  if(aB$intercept){
    A <- diag(1,k); A[,cB] <- -aB$M; A[cB, cB] <- 1
    A <- blockdiag(A,g,2)
    th <- A%*%th
    covar <- A%*%covar%*%t(A)
  }

  ##########################

  A <- blockdiag(aX$rot,g)
  th <- A%*%th
  covar <- A%*%covar%*%t(A)

  A <- blockdiag(t(aB$rot),g,2)
  th <- A%*%th
  covar <- A%*%covar%*%t(A)

  ##########################

  A <- blockdiag(mydiag(1/aX$s),g)
  th <- A%*%th
  covar <- A%*%covar%*%t(A)

  if(aX$intercept){
    A <- diag(1,q); A[cX,] <- -aX$m/aX$s[cX]; A[cX, cX] <- 1
    A <- blockdiag(A,g)
    th <- A%*%th
    covar <- A%*%covar%*%t(A)
  }

  ##########################

  A <- blockdiag(mydiag(1/aB$s),g,2)
  th <- A%*%th
  covar <- A%*%covar%*%t(A)

  if(aB$intercept){
    A <- diag(1,k); A[,cB] <- -aB$m/aB$s[cB]; A[cB, cB] <- 1
    A <- blockdiag(A,g,2)
    th <- A%*%th
    covar <- A%*%covar%*%t(A)
  }

  ##########################

  v <- (ay$M - ay$m)/10
  th <- th*v
  covar <- covar*(v^2)
  theta <- matrix(th,q,k)
  theta[cX,cB] <- theta[cX,cB] + ay$m/aB$s[cB]/aX$s[cX]

  list(theta = theta, covar = covar)
}

piqr.newton <- function(theta, y,z,d,X,Xw, bfun, s, type, tol=1e-6, maxit=100, safeit, eps0,
                       segno=1, lambda=0){

  p.bisec.internal <- getFromNamespace("p.bisec.internal", "qrcm")

  if(type == "iqr"){ee <- piqr.ee}
  else if(type == "ciqr"){ee <- pciqr.ee}
  else{ee <- pctiqr.ee}

  q <- nrow(theta)
  k <- ncol(theta)
  s <- c(s == 1)

  p.star.y <- p.bisec.internal(theta, y,X, bfun$bp)
  if(type == "ctiqr"){p.star.z <- p.bisec.internal(theta, z,X, bfun$bp)}
  G <- ee(theta, y,z,d,X,Xw, bfun, p.star.y, p.star.z, J = FALSE, segno=segno, lambda=lambda)

  g <- G$g[s]
  conv <- FALSE
  eps <- eps0

  for(i in 1:safeit){

    if(conv | max(abs(g)) < tol){break}
    u <- rep.int(0, q*k)
    u[s] <- g
    delta <- matrix(u, q,k)
    delta[is.na(delta)] <- 0
    cond <- FALSE

    while(!cond){

      new.theta <- theta - delta*eps
      if(max(abs(delta*eps)) < tol){conv <- TRUE; break}
      p.star.y <- p.bisec.internal(new.theta, y,X, bfun$bp)
      if(type == "ctiqr"){p.star.z <- p.bisec.internal(new.theta, z,X, bfun$bp)}
      G1 <- ee(new.theta, y,z,d,X,Xw, bfun, p.star.y, p.star.z, J = FALSE, segno=segno, lambda=lambda)
      g1 <- G1$g[s]
      cond <- (sum(g1^2) < sum(g^2))
      eps <- eps*0.5
    }

    if(conv){break}
    g <- g1
    G <- G1
    theta <- new.theta
    eps <- min(eps*2,0.1)
  }

  # Newton-Raphson

  alg <- "nr"
  conv <- FALSE
  eps <- 0.1
  h <- ee(theta, y,z,d,X,Xw, bfun, p.star.y, p.star.z, J = TRUE, G = G, segno=segno, lambda=lambda)$J[s,s, drop = FALSE]
  h <- h + diag(1e-4, nrow(h))

  for(i in 1:maxit){

    if(conv | max(abs(g)) < tol){break}

    ####

    if(type == "iqr"){
      H1 <- try(chol(h), silent = TRUE)
      err <- (class(H1) == "try-error")
    }
    else{
      H1 <- qr(h)
      r <- H1$rank
      err <- (r != ncol(h))
    }
    if(!err){
      if(alg == "gs"){alg <- "nr"; eps <- 1}
      delta <- (if(type == "iqr") chol2inv(H1)%*%g else qr.solve(H1)%*%g)
    }
    else{
      if(alg == "nr"){alg <- "gs"; eps <- 1}
      delta <- g
    }

    u <- rep.int(0, q*k)
    u[s] <- delta
    delta <- matrix(u, q,k)
    delta[is.na(delta)] <- 0
    cond <- FALSE
    while(!cond){
      new.theta <- theta - delta*eps
      if(max(abs(delta*eps)) < tol){conv <- TRUE; break}
      p.star.y <- p.bisec.internal(new.theta, y,X, bfun$bp)
      if(type == "ctiqr"){p.star.z <- p.bisec.internal(new.theta, z,X, bfun$bp)}
      G1 <- ee(new.theta, y,z,d,X,Xw, bfun, p.star.y, p.star.z, J = FALSE, segno=segno, lambda=lambda)
      g1 <- G1$g[s]
      cond <- (sum(g1^2) < sum(g^2))
      eps <- eps*0.5
    }

    if(conv){break}
    g <- g1
    G <- G1
    theta <- new.theta
    h <- ee(theta, y,z,d,X,Xw, bfun, p.star.y, p.star.z, J = TRUE, G = G, segno=segno, lambda=lambda)$J[s,s, drop = FALSE]
    h <- h + diag(1e-4, nrow(h))

    if(i > 1){eps <- min(eps*10,1)}
    else{eps <- min(eps*10,0.1)}
  }

  p.star.y <- p.bisec.internal(theta, y,X, bfun$bp)
  py <- bfun$p[p.star.y]
  if(type == "ctiqr"){
    p.star.z <- p.bisec.internal(theta, z,X, bfun$bp)
    pz <- bfun$p[p.star.z]
    pz <- pmin(pz, py - 1e-8)
    pz[p.star.z == 1] <- 0
  }
  else{p.star.z <- pz <- NULL}

  G <- ee(theta, y, z, d, X, Xw, bfun, p.star.y, p.star.z, J=FALSE, segno=1, lambda=0)

  list(coefficients = matrix(theta, q, k),
       converged = (i < maxit), n.it = i, p.star.y = p.star.y, p.star.z = p.star.z, py = py, pz = pz,
       ee = g, jacobian = h, rank = (alg == "nr")*sum(s), grad=G$g)
}

piqr.ee <- function(theta, y, z, d, X, Xw, bfun, p.star.y, p.star.z, J=TRUE, G,
                    i=FALSE, segno=1, lambda=0){

  k <- ncol(theta)
  n <- length(y)
  BB1 <- bfun$BB1
  if(missing(G)){
    B <- bfun$Bp[p.star.y,, drop = FALSE]
    S1 <- BB1 - B
    if(!i){g <- c(crossprod(Xw,S1) - segno*lambda)}
    else{g <- NULL; for(h in 1:k){g <- cbind(g,X*S1[,h])}
    if(length(segno) != 1) g <- g - matrix(rep(c(segno)*lambda/n, n), nrow=n, byrow=T)}
  }
  else{B <- G$B; g <- G$g}


  if(J){
    # b <- bfun$bp[p.star.y,, drop = FALSE]
    b1 <- bfun$b1p[p.star.y,, drop = FALSE]
    bij <- bfun$bpij[p.star.y,, drop = FALSE]

    A1 <- 1/c(.rowSums(tcrossprod(X, t(theta))*b1, n,k))
    A1 <- pmax0(A1)
    A1[attr(p.star.y, "out")] <- 0
    Xw <- Xw*A1

    J <- NULL
    count <- 0
    for(i1 in 1:k){
      h.temp <- NULL
      for(i2 in 1:k){
        count <- count + 1
        h.temp <- cbind(h.temp, crossprod(Xw, X*bij[,count]))
      }
      J <- rbind(J, h.temp)
    }
  }

  list(g = g, J = J, B = B)
}

pciqr.ee <- function(theta, y, z, d, X, Xw, bfun, p.star.y, p.star.z, J = TRUE, G,
                     i = FALSE, segno=1, lambda=0){

  k <- ncol(theta)
  n <- length(y)
  BB1 <- bfun$BB1

  if(missing(G)){
    B <- bfun$Bp[p.star.y,, drop = FALSE]
    BB <- bfun$BBp[p.star.y,, drop = FALSE]
    py <- bfun$p[p.star.y]

    a <- (1 - d)/(1 - py); a[attr(p.star.y, "out.r")] <- 0
    S1 <- a*(BB - BB1)
    a <- d; a[attr(p.star.y, "out.r")] <- 1
    S1 <- BB1 - a*B + S1

    if(!i){g <- c(crossprod(Xw,S1) - segno*lambda)}
    else{g <- NULL; for(h in 1:k){g <- cbind(g,X*S1[,h])}
    if(length(segno) != 1) g <- g - matrix(rep(c(segno)*lambda/n, n), nrow=n, byrow=T)}
  }
  else{B <- G$B; BB <- G$BB; g <- G$g; py <- G$py}

  if(J){
    b <- bfun$bp[p.star.y,, drop = FALSE]
    b1 <- bfun$b1p[p.star.y,, drop = FALSE]

    A1 <- 1/c(.rowSums(tcrossprod(X, t(theta))*b1, n,k))
    A1 <- pmax0(A1)
    A1[attr(p.star.y, "out")] <- 0
    a <- 1 - py
    a[attr(p.star.y, "out.r")] <- 1
    # the value 1 is arbitrary, only to avoid NAs (will be canceled by the zeroes in A1)
    A2 <- d*b + (1 - d)/(a^2)*(BB1 - B*a - BB)
    Xw <- Xw*A1

    J <- NULL
    for(i1 in 1:k){
      h.temp <- NULL
      for(i2 in 1:k){h.temp <- cbind(h.temp, crossprod(Xw, X*(b[,i2]*A2[,i1])))}
      J <- rbind(J, h.temp)
    }
  }
  list(g = g, J = J, B = B, BB = BB, py = py)
}

pctiqr.ee <- function(theta, y, z, d, X, Xw, bfun, p.star.y, p.star.z, J = TRUE, G,
                      i = FALSE, segno=1, lambda=0){

  k <- ncol(theta)
  n <- length(y)
  BB1 <- bfun$BB1
  if(missing(G)){
    B.y <- bfun$Bp[p.star.y,, drop = FALSE]
    BB.y <- bfun$BBp[p.star.y,, drop = FALSE]
    B.z <- bfun$Bp[p.star.z,, drop = FALSE]
    BB.z <- bfun$BBp[p.star.z,, drop = FALSE]
    py <- bfun$p[p.star.y]
    pz <- bfun$p[p.star.z]

    out.y <- attr(p.star.y, "out.r")
    out.z <- attr(p.star.z, "out.r")
    a <- d

    S1.y <- (1 - d)/(1 - py)*(BB.y - BB1)
    S1.z <- 1/(1 - pz)*(BB.z - BB1)
    S1.y[out.y,] <- 0; a[out.y] <- 1
    S1.y[out.z,] <- S1.z[out.z,] <- a[out.z] <- 0

    S1 <- -a*B.y + S1.y - S1.z
    if(!i){g <- c(crossprod(Xw,S1) - segno*lambda)}
    else{g <- NULL; for(h in 1:k){g <- cbind(g,X*S1[,h])}
    if(length(segno) != 1) g <- g - matrix(rep(c(segno)*lambda/n, n), nrow=n, byrow=T)}
  }
  else{B.y <- G$B.y; BB.y <- G$BB.y; B.z <- G$B.z; BB.z <- G$BB.z; g <- G$g; py <- G$py; pz <- G$pz}

  if(J){
    b.y <- bfun$bp[p.star.y,, drop = FALSE]
    b1.y <- bfun$b1p[p.star.y,, drop = FALSE]
    b.z <- bfun$bp[p.star.z,, drop = FALSE]
    b1.z <- bfun$b1p[p.star.z,, drop = FALSE]

    Xtheta <- tcrossprod(X, t(theta))

    A1.y <- 1/c(.rowSums(Xtheta*b1.y, n,k))
    A1.z <- 1/c(.rowSums(Xtheta*b1.z, n,k))
    A1.y <- pmax0(A1.y)
    A1.z <- pmax0(A1.z)
    A1.y[attr(p.star.y, "out")] <- 0
    A1.z[attr(p.star.z, "out")] <- 0

    ay <- 1 - py; az <- 1 - pz
    ay[attr(p.star.y, "out.r")] <- az[attr(p.star.z, "out.r")] <- 1
    # the value 1 is arbitrary, only to avoid NAs (will be canceled by the zeroes in A1.y and A1.z)

    A2.y <- d*b.y - (1 - d)/(ay^2)*(B.y*ay + BB.y - BB1)
    A2.z <- 1/(az^2)*(B.z*az + BB.z - BB1)
    H.y <- A2.y*A1.y
    H.z <- A2.z*A1.z

    J <- NULL
    for(i1 in 1:k){
      h.temp <- NULL
      for(i2 in 1:k){h.temp <- cbind(h.temp, crossprod(Xw, X*(b.y[,i2]*H.y[,i1] + b.z[,i2]*H.z[,i1])))}
      J <- rbind(J, h.temp)
    }
  }
  list(g = g, J = J, B.y = B.y, BB.y = BB.y, B.z = B.z, BB.z = BB.z, py = py, pz = pz)
}

cov.theta.2 <- function(theta, y,z,d,X,Xw, weights, bfun,
                      p.star.y, p.star.z, type, s, segno=1, lambda=0, J){

  if(type == "iqr"){ee <- piqr.ee}
  else if(type == "ciqr"){ee <- pciqr.ee}
  else{ee <- pctiqr.ee}
  s <- c(s == 1)

  G.i <- ee(theta, y,z,d,X,Xw, bfun, p.star.y, p.star.z, J = FALSE, i = TRUE, segno=segno, lambda=lambda)$g
  s.i <- G.i[,s, drop = FALSE]
  G <- J

  Omega <- chol2inv(chol(t(s.i*weights)%*%s.i))
  Q0 <- t(G)%*%Omega%*%G
  Q <- chol2inv(chol(Q0))
  U <- matrix(0, length(s), length(s))
  U[s,s] <- Q

  list(Q = U, Q0 = Q0, jacobian = G, ee = colSums(s.i*weights), Omega = Omega, s.i = s.i)
}

#' @export
print.piqr <- function(x, digits=max(3L, getOption("digits") - 3L), ...){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  print.default(cbind(Df = x$df, ObjFunction = signif(x$obj.function, digits),
                      Lambda = signif(x$lambda, digits)), print.gap = 2L, quote = FALSE)

  cat("\n")
  invisible(x)
}

#' @export
plot.piqr <- function(x, xvar=c("norm", "lambda", "objective", "grad", "beta"), label=FALSE, lambda, which=NULL,
                      ask=TRUE, polygon=TRUE, ...){
  xvar = match.arg(xvar)
  if(xvar == "beta" & missing(lambda)){
    stop("Please insert a lambda value")
  }else{
    if(xvar == "beta" & !missing(lambda)){
      l <- match(lambda, x$lambda)
      if(is.na(l)) stop("Lambda value is not correct")
    }
  }

  if(xvar != "beta"){
    lambda <- NULL
    lambda <- x$lambda
    nl <- x$internal$nlambda
    intercept <- x$internal$intercept
    k <- x$internal$k
    p <- ncol(x$internal$X)
    if(intercept){
      beta <- sapply(1:nl, function(.i) x$coefficients[[.i]][-1, ])
      grad <- sapply(1:nl, function(.i) c(matrix(x$dl[, .i], ncol=k)[-1, ]))
    }else{
      beta <- sapply(1:nl, function(.i) x$coefficients[[.i]])
    }

    plotCoef(beta, lambda=lambda, df=x$df, dev=x$obj.function, grad=x$dl, label=label, xvar=xvar, ...)
  }else{
    plot.piqr2(x, lambda, which=which, ask=ask, polygon=polygon, ...)
  }
}

plot.piqr2 <- function(x, lambda, conf.int=TRUE, polygon=TRUE, which=NULL, ask=TRUE, ...){

  plot.piqr.int <- function(p,u,j,conf.int,L){
    beta <- u[[j]]$beta
    if(is.null(L$ylim)){
      if(conf.int){y1 <- min(u[[j]]$low); y2 <- max(u[[j]]$up)}
      else{y1 <- min(beta); y2 <- max(beta)}
      L$ylim <- c(y1,y2)
    }
    plot(p, u[[j]]$beta, xlab = L$xlab, ylab = L$ylab[j], main = L$labels[j],
         type = "l", lwd = L$lwd, xlim = L$xlim, ylim = L$ylim, col = L$col, cex.lab = L$cex.lab, cex.axis = L$cex.axis)
    if(conf.int){
      if(polygon){
        yy <- c(u[[j]]$low, tail(u[[j]]$up, 1), rev(u[[j]]$up), u[[j]]$low[1])
        xx <- c(p, tail(p, 1), rev(p), p[1])
        polygon(xx, yy, col = adjustcolor(L$col, alpha.f = 0.25), border = NA)
      }else{
        points(p, u[[j]]$low, lty = 2, lwd = L$lwd, type = "l", col = L$col, cex.lab = L$cex.lab, cex.axis = L$cex.axis)
        points(p, u[[j]]$up, lty = 2, lwd = L$lwd, type = "l", col = L$col, cex.lab = L$cex.lab, cex.axis = L$cex.axis)
      }
    }
    # if(conf.int){
    #   points(p, u[[j]]$low, lty = 2, lwd = L$lwd, type = "l", col = L$col, cex.lab = L$cex.lab, cex.axis = L$cex.axis)
    #   points(p, u[[j]]$up, lty = 2, lwd = L$lwd, type = "l", col = L$col, cex.lab = L$cex.lab, cex.axis = L$cex.axis)
    # }
  }

  L <- list(...)
  if(is.null(L$xlim)){L$xlim = c(0.01,0.99)}
  if(is.null(L$lwd)){L$lwd <- 2}
  if(is.null(L$cex.lab)){L$cex.lab <- 1}
  if(is.null(L$cex.axis)){L$cex.axis <- 1}
  if(is.null(L$col)){L$col <- "black"}
  if(is.null(L$xlab)){L$xlab <- "p"}
  if(is.null(L$labels)) L$labels <- rownames(x$coefficients[[length(x$lambda)]])
  q <- length(L$labels)
  if(is.null(L$ylab)){L$ylab <- rep("beta(p)", q)}
  L$labels <- c(L$labels, "qqplot")

  p <- seq.int(max(0.001,L$xlim[1]), min(0.999,L$xlim[2]), length.out = 100)
  u <- predict.piqr(x, lambda, p=p, type="beta", se=conf.int)

  if(!is.null(which) | !ask){
    if(is.null(which)){which <- 1:q}
    for(j in which){plot.piqr.int(p,u,j,conf.int,L)}
  }
  else{
    pick <- 1
    while(pick > 0 && pick <= q + 1){
      pick <- menu(L$labels, title = "Make a plot selection (or 0 to exit):\n")
      if(pick > 0 && pick <= q){plot.piqr.int(p,u,pick,conf.int,L)}
      else if(pick == q + 1){
        KM <- NULL
        KM$time <- sort(runif(nrow(x$internal$mf), 0, 1))
        KM$cdf <- sort(predict(x, lambda, type="CDF")[, 1])
        plot(KM$time, KM$cdf, pch = 20, cex = 0.5,
             xlim = c(0,1), ylim = c(0,1),
             ylab = "U(0,1) quantiles", xlab = "fitted CDF quantiles", cex.axis = L$cex.axis, cex.lab = L$cex.lab)
        abline(0,1)
      }
    }
  }
}

#' @export
predict.piqr <- function(object, lambda, type = c("beta", "CDF", "QF", "sim"), newdata, p, se=TRUE, ...){
  if(missing(lambda)) stop("Please insert a lambda value")
  l <- match(lambda, object$lambda)
  if(is.na(l)) stop("Lambda value is not correct")

  pred.beta <- getFromNamespace("pred.beta", "qrcm")
  p.bisec <- getFromNamespace("p.bisec", "qrcm")
  apply_bfun <- getFromNamespace("apply_bfun", "qrcm")
  extract.p <- getFromNamespace("extract.p", "qrcm")

  if(is.na(match(type <- type[1], c("beta", "CDF", "QF", "sim")))){stop("invalid 'type'")}
  if(type == "beta"){
    if(missing(p)){p <- seq.int(0.01,0.99,0.01)}
    if(any(p <= 0 | p >= 1)){stop("0 < p < 1 is required")}
    model <- list(coefficients=object$coefficient[[l]], mf=object$internal$mf, covar=object$covar[[l]])
    return(pred.beta(model, p, se=se))
  }

  mf <- object$internal$mf
  mt <- terms(mf)
  miss <- attr(mf, "na.action")
  fittype <- attr(mf, "type")
  nomiss <- (if(is.null(miss)) 1:nrow(mf) else (1:(nrow(mf) + length(miss)))[-miss])
  xlev <- .getXlevels(mt, mf)

  if(!missing(newdata)){
    if(type == "CDF"){
      yn <- as.character(if(fittype == "ctiqr") mt[[2]][[3]]
                         else if(fittype == "ciqr") mt[[2]][[2]] else mt[[2]])
      if(is.na(ind <- match(yn, colnames(newdata))))
      {stop("for 'type = CDF', 'newdata' must contain the y-variable")}
      if(fittype == "ciqr"){newdata[,as.character(mt[[2]][[3]])] <- 1}
      if(fittype == "ctiqr"){newdata[,as.character(mt[[2]][[4]])] <- 1
      newdata[,as.character(mt[[2]][[2]])] <- -Inf}
    }else{mt <- delete.response(mt)}
    if(any(is.na(match(all.vars(mt), colnames(newdata)))))
    {stop("'newdata' must contain all x-variables")}

    mf <- model.frame(mt, data = newdata, xlev = xlev)
    if(nrow(mf) == 0){
      nr <- nrow(newdata)
      if(type == "CDF"){
        out <- data.frame(matrix(NA,nr,3))
        colnames(out) <- c("log.f", "log.F", "log.S")
        rownames(out) <- rownames(newdata)
      }else if(type == "QF"){
        out <- data.frame(matrix(NA,nr,length(p)))
        colnames(out) <- paste("p",p, sep = "")
        rownames(out) <- rownames(newdata)
      }else{out <- rep.int(NA, nr)}
      return(out)
    }
    miss <- attr(mf, "na.action")
    nomiss <- (if(is.null(miss)) 1:nrow(mf) else (1:nrow(newdata))[-miss])
  }

  x <- model.matrix(mt, mf)

  if(type == "CDF"){
    bfun <- attr(object$internal$mf, "bfun")
    y <- cbind(model.response(mf))[,1 + (fittype == "ctiqr")]
    Fy <- p.bisec(object$coefficients[[l]], y,x, bfun)
    b1 <- apply_bfun(bfun, Fy, "b1fun")
    fy <- 1/c(rowSums((x%*%object$coefficients[[l]])*b1))
    fy[attr(Fy, "out")] <- 0
    if(any(fy < 0)){warning("some PDF values are negative (quantile crossing)")}
    CDF <- PDF <- NULL
    CDF[nomiss] <- Fy
    PDF[nomiss] <- fy
    CDF[miss] <- PDF[miss] <- NA
    out <- data.frame(CDF = CDF, PDF = PDF)
    rownames(out)[nomiss] <- rownames(mf)
    if(!is.null(miss)){rownames(out)[miss] <- names(miss)}
    return(out)
  }
  else if(type == "QF"){
    if(missing(p)){stop("please indicate the value(s) of 'p' to compute x*beta(p)")}
    if(any(p <= 0 | p >= 1)){stop("0 < p < 1 is required")}

    fit <- se.fit <- matrix(NA, length(c(miss,nomiss)), length(p))
    model <- list(coefficients=object$coefficients[[l]], covar=object$covar[[l]], mf=object$internal$mf)
    for(j in 1:length(p)){
      fit.beta <- extract.p(model, p[j], cov = se)
      fit[nomiss,j] <- x%*%cbind(fit.beta$coef[,1])
      if(se){se.fit[nomiss,j] <- sqrt(diag(x%*%fit.beta$cov%*%t(x)))}
    }
    fit <- data.frame(fit)
    colnames(fit) <- paste("p",p, sep = "")
    rownames(fit)[nomiss] <- rownames(mf)
    if(!is.null(miss)){rownames(fit)[miss] <- names(miss)}
    if(se){
      se.fit <- data.frame(se.fit)
      colnames(se.fit) <- paste("p",p, sep = "")
      rownames(se.fit)[nomiss] <- rownames(mf)
      if(!is.null(miss)){rownames(se.fit)[miss] <- names(miss)}
      return(list(fit = fit, se.fit = se.fit))
    }
    else{return(fit)}
  }
  else{
    p <- runif(nrow(x))
    beta <- apply_bfun(attr(object$internal$mf, "bfun"), p, "bfun")%*%t(object$coefficients[[l]])
    y <- NULL; y[nomiss] <- rowSums(beta*x); y[miss] <- NA
    return(y)
  }
}

#' @export
summary.piqr <- function(object, lambda, SE=FALSE, p, cov=FALSE, ...){
  if(missing(lambda)) stop("Please insert a lambda value")
  l <- match(lambda, object$lambda)
  if(is.na(l)) stop("Lambda value is not correct")

  if(missing(p)){
    covar <- object$covar[[l]]
    if(SE){
      se <- sqrt(diag(covar))
      se <- matrix(se, ncol=object$internal$k,
                   dimnames=list(attr(object$internal$mf, "stats")$X$coef.names,
                                 attr(object$internal$mf, "stats")$B$coef.names))
    }else{
      se <- NULL
    }

    out <- list(nlambda=lambda, coefficients=object$coefficients[[l]], se=se, covar=covar,
                obj.function=object$obj.function[l], n=nrow(object$internal$mf), free.par=object$df[l],
                call=object$call, SE=SE)
  }else{
    extract.p <- getFromNamespace("extract.p", "qrcm")
    apply_bfun <- getFromNamespace("apply_bfun", "qrcm")

    model <- list(coefficients=object$coefficients[[l]], covar=object$covar[[l]], mf=object$internal$mf)

    out <- list()
    for(i in 1:length(p)){
      out[[i]] <- extract.p(model, p[i], cov)
    }
    names(out) <- paste("p =", p)
    attr(out, "nacoef") <- which(apply(model$coefficients, 1, function(v){all(v == 0)}))
    out$call <- object$call
  }
  class(out) <- "summary.piqr"

  return(out)
}

#' @export
print.summary.piqr <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  cat("\nCall: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if(!is.null(x$coef)){
    nacoef <- which(x$coef == 0)
    x$coef[nacoef] <- NA
    x$se[nacoef] <- NA

    cat("Selected value of lambda:", format(x$nlambda, digits = digits), "\n")
    cat("n. of observations:", x$n, "\n")
    cat("n. of free parameters:", x$free.par, "\n\n")

    cat("######################", "\n")
    cat("######################", "\n\n")

    cat("Coefficients:\n")
    print.default(format(x$coef, digits = digits), print.gap = 2L, quote = FALSE)
    cat("\n")
    if(x$SE){
      cat("Standard errors:\n")
      print.default(format(x$se, digits = digits), print.gap = 2L,
                    quote = FALSE)
      cat("\n")
    }

    cat("######################", "\n")
    cat("######################", "\n\n")

    cat("Minimized loss function:", x$obj.function)
    cat("\n\n")
  }else{
    nacoef <- attr(x, "nacoef")
    for(j in 1:(length(x) - 1)){
      cat(paste(names(x)[j], "\n"))
      cat("\n")
      cat("Coefficients:\n")
      coe <- x[[j]]$coef; coe[nacoef,] <- NA
      printCoefmat(coe, digits = digits, signif.stars = TRUE, cs.ind = 1:2, tst.ind = 3,
                   P.values = TRUE, has.Pvalue = TRUE)
      cat("\n")

      if(!is.null(x[[j]]$cov)){
        cat("Covar:\n")
        print.default(format(x[[j]]$cov, digits = digits), print.gap = 2L, quote = FALSE)
      }
      cat("\n\n")
    }
  }

  invisible(x)
}

plotCoef <- function (beta, norm, lambda, df, dev, grad, label=FALSE,
                      xvar=c("norm", "lambda", "objective", "grad"), xlab=iname, ylab="Coefficients", ...){
  which = nonzeroCoef(beta)
  nwhich = length(which)
  switch(nwhich + 1, `0` = {
    warning("No plot produced since all coefficients zero")
    return()
  }, `1` = warning("1 or less nonzero coefficients; glmnet plot is not meaningful"))
  beta = as.matrix(beta[which, , drop = FALSE])
  xvar = match.arg(xvar)
  switch(xvar, norm = {
    index = if (missing(norm)) apply(abs(beta), 2, sum) else norm
    iname = "L1 Norm"
    approx.f = 1
  }, lambda = {
    index = log(lambda)
    iname = expression(log(lambda))
    approx.f = 0
  }, objective = {
    index = dev
    iname = "Objective Function"
    approx.f = 1
  }, grad = {
    index = log(lambda)
    iname = expression(log(lambda))
    beta = grad
    ylab = "Gradient"
    approx.f = 0
  })
  dotlist = list(...)
  type = dotlist$type
  if(is.null(type))
    matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab, type = "l", ...)
  else matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab, ...)
  atdf = pretty(index)
  prettydf = approx(x = index, y = df, xout = atdf, rule = 2,
                    method = "constant", f = approx.f)$y
  axis(3, at = atdf, labels = prettydf, tcl = NA, ...)
  if (label) {
    nnz = length(which)
    xpos = max(index)
    pos = 4
    if (xvar == "lambda") {
      xpos = min(index)
      pos = 2
    }
    xpos = rep(xpos, nnz)
    ypos = beta[, ncol(beta)]
    text(xpos, ypos, paste(which), cex = 0.5, pos = pos)
  }
}

nonzeroCoef <- function (beta, bystep=FALSE){
  nr = nrow(beta)
  if (nr == 1) {
    if (bystep)
      apply(beta, 2, function(x) if (abs(x) > 0)
        1
        else NULL)
    else {
      if (any(abs(beta) > 0))
        1
      else NULL
    }
  }
  else {
    beta = abs(beta) > 0
    which = seq(nr)
    ones = rep(1, ncol(beta))
    nz = as.vector((beta %*% ones) > 0)
    which = which[nz]
    if (bystep) {
      if (length(which) > 0) {
        beta = as.matrix(beta[which, , drop = FALSE])
        nzel = function(x, which) if (any(x))
          which[x]
        else NULL
        which = apply(beta, 2, nzel, which)
        if (!is.list(which))
          which = data.frame(which)
        which
      }
      else {
        dn = dimnames(beta)[[2]]
        which = vector("list", length(dn))
        names(which) = dn
        which
      }
    }
    else which
  }
}

