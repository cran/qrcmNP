\name{summary.niqr}
\alias{summary.niqr}
\title{
Summary After Nonlinear Quantile Regression Coefficients Modeling
}
\description{
Summary of an object of class \dQuote{\code{niqr}}.
}
\usage{
\method{summary}{niqr}(object, p, \ldots)
}
\arguments{
  \item{object}{
  an object of class \dQuote{\code{niqr}}, the result of a call to \code{\link{niqr}}.
}
  \item{p}{
  an optional vector of quantiles.
}
  \item{\ldots}{
  for future methods.
}
}
\details{
A summary of the model is printed.
}
\author{
Gianluca Sottile \email{gianluca.sottile@unipa.it}
}
\seealso{
\code{\link{niqr}}, for model fitting; \code{\link{testfit.niqr}}, for goodness of fit test; \code{\link{predict.niqr}} and \code{\link{plot.niqr}}, for predicting and plotting objects of class \dQuote{\code{niqr}}.
}
\examples{

n <- 300
x <- runif(n)
fun <- function(theta, p){
  beta0 <- theta[1] + exp(theta[2]*p)
  beta1 <- theta[3] + theta[4]*p
  cbind(beta0, beta1)}
beta <- fun(c(1,1,1,1), runif(n))
y <- beta[, 1] + beta[, 2]*x
model <- niqr(fun=fun, x0=rep(0, 4), X=cbind(1,x), y=y)

summary(model)
summary(model, p=c(.01,.05))
}
