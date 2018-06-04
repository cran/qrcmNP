\name{plot.piqr}
\alias{plot.piqr}
\title{
Plot Penalized Quantile Regression Coefficients
}
\description{
Produces a coefficient profile plot of the quantile
regression coefficient paths for a fitted model of
class \dQuote{\code{piqr}}.
}
\usage{
\method{plot}{piqr}(x, xvar=c("norm", "lambda", "objective", "grad", "beta"), label=FALSE,
                    lambda, which=NULL, ask=TRUE, polygon=TRUE, ...)
}
\arguments{
  \item{x}{
  an object of class \dQuote{\code{piqr}}, typically the result of a call to \code{\link{piqr}}.
}
  \item{xvar}{
  What is on the X-axis. "norm" plots against the L1-norm of the coefficients,
  "lambda" against the log-lambda sequence, "objective" against the value
  of the minimized integrated loss function and "grad" the log-lambda sequence
  against the gradient.
  xvar = "beta" needs a lambda value to plot quantile regression coefficients
  \eqn{\beta(p | \theta(\lambda))} as a function of p, based on the fitted model of class \dQuote{\code{piqr}}
}
  \item{label}{
  If TRUE, label the curves with variable sequence numbers.
}
  \item{lambda}{
  a value of lambda in the sequence of the object of class \dQuote{\code{piqr}}. Could be the best
  after selecting the result of a call to \code{\link{gof.piqr}}
}
  \item{which}{
  an optional numerical vector indicating which coefficient(s) to plot. If which = NULL, all coefficients are     plotted.
}
  \item{ask}{
  logical. If which = NULL and ask = TRUE (the default), you will be asked interactively which coefficients to plot.
}
\item{polygon}{
  ogical. If TRUE, confidence intervals are represented by shaded areas via polygon. Otherwise, dashed lines are used.
}
  \item{\ldots}{
  additional graphical parameters, that can include \kbd{xlim, ylim, xlab, ylab, col, lwd}.
  See \code{\link{par}}.}
}
\details{
A coefficient profile plot is produced.
}
\author{
Gianluca Sottile \email{gianluca.sottile@unipa.ot}
}
\seealso{
\code{\link{piqr}} for model fitting; \code{\link{gof.piqr}} for the model selection criteria; \code{\link{summary.piqr}} and \code{\link{predict.piqr}} for model summary and prediction.
}
\examples{

  # using simulated data

  n <- 300
  x <- runif(n)
  qy <- function(p,x){p^2 + x*log(p)}
  # true quantile function: Q(p | x) = beta0(p) + beta1(p)*x, with
     # beta0(p) = p^2
     # beta1(p) = log(p)
  y <- qy(runif(n), x) # to generate y, plug uniform p in qy(p,x)

  obj <- piqr(y ~ x, formula.p = ~ slp(p,3), nlambda=50)
  best <- gof.piqr(obj, method="BIC", plot=FALSE)
  par(mfrow = c(2,2))
  plot(obj, xvar="norm")
  plot(obj, xvar="lambda")
  plot(obj, xvar="objective")
  plot(obj, xvar="grad")
  par(mfrow=c(1,2));plot(obj, xvar="beta", lambda=best$minLambda, ask=FALSE)
  # flexible fit with shifted Legendre polynomials
}
\keyword{methods}

