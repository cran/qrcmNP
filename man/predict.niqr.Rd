\name{predict.niqr}
\alias{predict.niqr}
\title{
Prediction After Nonlinear Quantile Regression Coefficients Modeling
}
\description{
Predictions from an object of class \dQuote{\code{niqr}}.
}
\usage{
\method{predict}{niqr}(object, type=c("beta", "CDF", "QF", "sim"), newdata, p, \ldots)
}
\arguments{
  \item{object}{
  an object of class \dQuote{\code{niqr}}, the result of a call to \code{\link{niqr}}.
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
\item{\ldots}{
  for future methods.
}
}
\details{
Different type of prediction from the model.
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
\code{\link{niqr}}, for model fitting; \code{\link{test.fit.niqr}}, to do goodness of fit test; \code{\link{summary.niqr}} and \code{\link{plot.niqr}}, for summarizing and plotting \code{niqr} objects.
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
model <- niqr(fun=fun, x0=rep(0, 4), X=cbind(1,x), y=y)

# predict beta(0.25), beta(0.5), beta(0.75)
predict(model, type = "beta", p = c(0.25,0.5, 0.75))

# predict the CDF and the PDF at new values of x and y
predict(model, type = "CDF",
        newdata = data.frame(X1=runif(3), y = c(1,2,3)))

# computes the quantile function at new x, for p = (0.25,0.5,0.75)
predict(model, type = "QF", p = c(0.25,0.5,0.75),
        newdata = data.frame(X1=runif(3), y = c(1,2,3)))

# simulate data from the fitted model
ysim <- predict(model, type = "sim") # 'newdata' can be supplied

# if the model is correct, the distribution of y and that of ysim should be similar
qy <- quantile(y, prob = seq(.1,.9,.1))
qsim <- quantile(ysim, prob = seq(.1,.9,.1))
plot(qy, qsim); abline(0,1)
}
