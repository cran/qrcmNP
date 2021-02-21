\name{summary.piqr}
\alias{summary.piqr}
\title{
Summary After Penalized Quantile Regression Coefficients Modeling
}
\description{
Summary of an object of class \dQuote{\code{piqr}}, after selecting the best tuning parameter.
}
\usage{
\method{summary}{piqr}(object, pos.lambda, SE=FALSE, p, cov=FALSE, \ldots)
}
\arguments{
  \item{object}{
  an object of class \dQuote{\code{piqr}}, the result of a call to \code{\link{piqr}}.
}
  \item{pos.lambda}{
  the position of a lambda in the sequence of the object of class \dQuote{\code{piqr}}. Could be the best
  after selecting the result of a call to \code{\link{gof.piqr}}
}
  \item{SE}{
  if TRUE standard errors are printed. Standard errors are computed through sandwich formula only for the             regularized parameters.
}
  \item{p}{
  an optional vector of quantiles.
}
  \item{cov}{
  ff TRUE, the covariance matrix of \eqn{\beta(p)} is reported. Ignored if p is missing.
}
\item{\ldots}{for future methods.}
}
\details{
If the best lambda or one value of lambda is chosen a summary of the selected model is printed.
}
\value{
See details in \code{\link{summary.iqr}}
}
\author{
Gianluca Sottile \email{gianluca.sottile@unipa.it}
}
\seealso{
\code{\link{piqr}}, for model fitting; \code{\link{gof.piqr}}, to find the best lambda value; \code{\link{predict.piqr}} and \code{\link{plot.piqr}}, for predicting and plotting objects of class \dQuote{\code{piqr}}.
}
\examples{

# using simulated data

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
