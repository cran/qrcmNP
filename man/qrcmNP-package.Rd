\name{qrcmNP-package}
\alias{qrcmNP-package}
\docType{package}
\title{
Nonlinear and Penalized Quantile Regression Coefficients Modeling
}
\description{
This package implements a nonlinear Frumento and Bottai's (2016) method for quantile regression coefficient modeling (qrcm), in which quantile regression coefficients are described by (flexible) parametric functions of the order of the quantile. In the classical qrcm framework the linearity in \eqn{\beta(p)} and/or in \eqn{\theta} could be relaxed at a cost of more complicated expressions for the ojective and the gradient functions. Here, we propose an efficiently algorithm to use more flexible structures for the regression coefficients. With respect to the most famous function nlrq (quantreg package) our main function niqr implements the integrated quantile regression idea  of Frumento and Bottai's (2016) for nonlinear functions. As already known, this practice allows to estimate quantiles all at one time and not one at a time.
This package also implements a penalized Frumento and Bottai's (2015) method for quantile regression coefficient modeling (qrcm). This package fits lasso qrcm using pathwise coordinate descent algorithm. With respect to some other packages which implements the L1-quantile regression (e.g. quantreg, rqPen) estimating quantiles one at a time our proposal allows to estimate the conditional quantile function parametrically estimating quantiles all at one and to do variable selction in the meanwhile.
}
\details{
\tabular{ll}{
Package: \tab qrcmNP\cr
Type: \tab Package\cr
Version: \tab 0.2.1\cr
Date: \tab 2024-01-22\cr
License: \tab GPL-2\cr
}
The function \code{\link{niqr}} permits specifying nonlinear basis for each variables. The function \code{\link{testfit.niqr}} permits to do goodness of fit.
The auxiliary functions \code{\link{summary.niqr}}, \code{\link{predict.niqr}}, and \code{\link{plot.niqr}} can be used to extract information from the fitted model.
The function \code{\link{piqr}} permits specifying the lasso regression model. The function \code{\link{gof.piqr}} permits to select the best tuning parameter through AIC, BIC, GIC and GCV criteria. The auxiliary functions \code{\link{summary.piqr}}, \code{\link{predict.piqr}}, and \code{\link{plot.piqr}} can be used to extract information from the fitted model.
}
\author{
Gianluca Sottile

Maintainer: Gianluca Sottile <gianluca.sottile@unipa.it>
}
\keyword{ package }
\references{
Sottile G, Frumento P, Chiodi M, Bottai M. (2020). \emph{A penalized approach to covariate selection through quantile regression coefficient models}. Statistical Modelling, 20(4), pp 369-385. doi:10.1177/1471082X19825523.

Frumento, P., and Bottai, M. (2016). \emph{Parametric modeling of quantile regression coefficient functions}. Biometrics, 72(1), pp 74-84, doi:10.1111/biom.12410.

Friedman, J., Hastie, T. and Tibshirani, R. (2008). \emph{Regularization Paths for Generalized Linear Models via Coordinate Descent}. Journal of Statistical Software, Vol. 33(1), pp 1-22 Feb 2010.
}
\examples{

# use simulated data

n <- 300
x <- runif(n)
fun <- function(theta, p){
  beta0 <- theta[1] + exp(theta[2]*p)
  beta1 <- theta[3] + theta[4]*p
  cbind(beta0, beta1)}
beta <- fun(c(1,1,1,1), runif(n))
y <- beta[, 1] + beta[, 2]*x
model <- niqr(fun=fun, x0=rep(0, 4), X=cbind(1,x), y=y)

# use simulated data

set.seed(1234)
n <- 300
x1 <- rexp(n)
x2 <- runif(n, 0, 5)
x <- cbind(x1,x2)

b <- function(p){matrix(cbind(1, qnorm(p), slp(p, 2)), nrow=4, byrow=TRUE)}
theta <- matrix(0, nrow=3, ncol=4); theta[, 1] <- 1; theta[1,2] <- 1; theta[2:3,3] <- 2
qy <- function(p, theta, b, x){rowSums(x * t(theta \%*\% b(p)))}

y <- qy(runif(n), theta, b, cbind(1, x))

s <- matrix(1, nrow=3, ncol=4); s[1,3:4] <- 0
obj <- piqr(y ~ x1 + x2, formula.p = ~ I(qnorm(p)) + slp(p, 2), s=s, nlambda=50)

best <- gof.piqr(obj, method="AIC", plot=FALSE)
best2 <- gof.piqr(obj, method="BIC", plot=FALSE)

summary(obj, best$posMinLambda)
summary(obj, best2$posMinLambda)

}
