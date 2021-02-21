\name{plot.niqr}
\alias{plot.niqr}
\title{
Plot Nonlinear Quantile Regression Coefficients
}
\description{
Plots quantile regression coefficients \eqn{\beta(\theta, p)} as a function of p, based on a fitted model
of class \dQuote{\code{niqr}}.
}
\usage{
\method{plot}{niqr}(x, conf.int=TRUE, which=NULL, ask=TRUE, \ldots)
}
\arguments{
  \item{x}{
  an object of class \dQuote{\code{niqr}}, typically the result of a call to \code{\link{niqr}}.
}
  \item{conf.int}{
  logical. If TRUE, asymptotic 95\% confidence intervals are added to the plot.
}
  \item{which}{
  an optional numerical vector indicating which coefficient(s) to plot. If which = NULL,
  all coefficients are plotted.
}
\item{ask}{
  logical. If which = NULL and ask = TRUE (the default), you will be asked interactively
  which coefficients to plot.
}
  \item{\ldots}{
  additional graphical parameters, that can include \kbd{xlim, ylim, xlab, ylab, col, lwd}.
  See \code{\link{par}}.}
}
\author{
Gianluca Sottile \email{gianluca.sottile@unipa.ot}
}
\seealso{
\code{\link{niqr}} for model fitting; \code{\link{testfit.niqr}} for goodness of fit test; \code{\link{summary.niqr}} and \code{\link{predict.niqr}} for model summary and prediction.
}
\examples{

  # using simulated data

n <- 300
x <- runif(n)
fun <- function(theta, p){
  beta0 <- theta[1] + exp(theta[2]*p)
  beta1 <- theta[3] + theta[4]*p
  cbind(beta0, beta1)}
beta <- fun(c(1,1,1,1), runif(n))
y <- beta[, 1] + beta[, 2]*x
model <- niqr(fun=fun, x0=rep(0, 4), X=cbind(1, x), y=y)

plot(model, ask=FALSE)

}
\keyword{methods}

