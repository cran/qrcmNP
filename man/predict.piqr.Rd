\name{predict.piqr}
\alias{predict.piqr}
\title{
Prediction After Penalized Quantile Regression Coefficients Modeling
}
\description{
Predictions from an object of class \dQuote{\code{piqr}}, after selecting the best tuning parameter.
}
\usage{
\method{predict}{piqr}(object, lambda, type=c("beta", "CDF", "QF", "sim"), newdata,
                       p, se=TRUE, \ldots)
}
\arguments{
  \item{object}{
  an object of class \dQuote{\code{piqr}}, the result of a call to \code{\link{piqr}}.
}
  \item{lambda}{
  a value of lambda in the sequence of the object of class \dQuote{\code{piqr}}. Could be the best
  after selecting the result of a call to \code{\link{gof.piqr}}
}
  \item{type}{
  a character string specifying the type of prediction. See \sQuote{Details}.
}
  \item{newdata}{
  an optional data frame in which to look for variables with which to predict.
  If omitted, the data are used. For \kbd{type = "CDF"}, it must include the response variable.
 Ignored if \kbd{type = "beta"}.
}
  \item{p}{
  a numeric vector indicating the order(s) of the quantile to predict. Only used if
	\kbd{type = "beta"} or \kbd{type = "QF"}.
}
  \item{se}{
  logical. If TRUE (the default), standard errors of the prediction will be computed. Only used if type = "beta" or   type = "QF".
}
\item{\ldots}{for future methods.}
}
\details{
If the best lambda or one value of lambda is chosen, the function call \code{\link{predict.iqr}}.
}

\value{
See details in \code{\link{predict.iqr}}
}
\author{
Gianluca Sottile \email{gianluca.sottile@unipa.it}
}
\note{
Prediction may generate quantile crossing
if the support of the new covariates values supplied in \code{newdata}
is different from that of the observed data.
}
\seealso{
\code{\link{piqr}}, for model fitting; \code{\link{gof.piqr}}, to find the best lambda value; \code{\link{summary.piqr}} and \code{\link{plot.piqr}}, for summarizing and plotting \code{piqr} objects.
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

# predict beta(0.25), beta(0.5), beta(0.75)
predict(obj, best$minLambda, type = "beta", p = c(0.25,0.5, 0.75))

# predict the CDF and the PDF at new values of x and y
predict(obj, best$minLambda, type = "CDF",
        newdata = data.frame(x1=rexp(3), x2=runif(3), y = c(1,2,3)))

# computes the quantile function at new x, for p = (0.25,0.5,0.75)
predict(obj, best$minLambda, type = "QF", p = c(0.25,0.5,0.75),
        newdata = data.frame(x1=rexp(3), x2=runif(3), y = c(1,2,3)))

# simulate data from the fitted model
ysim <- predict(obj, best$minLambda, type = "sim") # 'newdata' can be supplied

# if the model is correct, the distribution of y and that of ysim should be similar
qy <- quantile(y, prob = seq(.1,.9,.1))
qsim <- quantile(ysim, prob = seq(.1,.9,.1))
plot(qy, qsim); abline(0,1)
}
