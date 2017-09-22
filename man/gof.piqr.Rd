\name{gof.piqr}
\alias{gof.piqr}
\title{
Goodness of Fit of Penalized Quantile Regression Coefficients Modeling
}
\description{
Goodness of Fit of an object of class \dQuote{\code{piqr}}, usefull to select the best tuning parameter.
}
\usage{
gof.piqr(object, method=c("AIC","BIC","GIC","GCV","NIC"), plot=TRUE, ...)
}
\arguments{
  \item{object}{
  an object of class \dQuote{\code{piqr}}, the result of a call to \code{\link{piqr}}.
}
  \item{method}{
  a method to evaluate the goodness of fit and select the best value of the tuning parameter.
}
  \item{plot}{
  if TRUE the chosen method is plotted - default is TRUE.
}
  \item{...}{
  additional arguments.
}
}
\details{
The best value of lambda is chosen minimizing the criterion, i.e. AIC, BIC, GIC, GCV.
}

\value{
\item{minLambda}{the best value of lambda.}
\item{dfMinLambda}{the number of nonzero parameters associated to the best lambda.}
\item{betaMin}{the parameters associated to the best lambda.}
\item{posMinLambda}{the position of the best lambda along the sequence of lambda.}
\item{call}{the matched call.}
}
\author{
Gianluca Sottile \email{gianluca.sottile@unipa.it}
}

\seealso{
\code{\link{piqr}}, for model fitting; \code{\link{summary.piqr}} and \code{\link{plot.piqr}}, for summarizing and plotting \code{piqr} objects.
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
obj <- piqr(y ~ x1 + x2, formula.p = ~ I(qnorm(p)) + slp(p, 2), s=s)

par(mfrow=c(1,2))
best <- gof.piqr(obj, method="AIC", plot=TRUE)
best2 <- gof.piqr(obj, method="BIC", plot=TRUE)
}
