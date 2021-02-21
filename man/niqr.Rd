\name{niqr}
\alias{niqr}
\title{
Nonlinear Quantile Regression Coefficients Modeling
}
\description{
This package implements a nonlinear Frumento and Bottai's (2015) method for quantile regression
coefficient modeling (qrcm), in which quantile regression coefficients are described
by (flexible) parametric functions of the order of the quantile.
}
\usage{
niqr(fun, fun2, x0, X, y, control=list())
}
\arguments{
  \item{fun}{
    a function of theta and p describing the beta functions.
  }
  \item{fun2}{
    a function of beta and X describing the whole quantile process.
  }
  \item{x0}{
    starting values to search minimum.
  }
  \item{X}{
    a design matrix containing the intercept.
  }
  \item{y}{
    the response variable.
  }
  \item{control}{
    a list of control parameters. See 'Details'.
  }
}
\details{
  Quantile regression permits modeling conditional quantiles of a response variabile,
  given a set of covariates.

  Assume that each coefficient can be expressed as a parametric function of \eqn{\theta, p} of the form:
  \deqn{\beta(\theta, p) = b_0(\theta_0, p) + b_1(\theta_1, p) + b_2(\theta_2 , p) + \ldots.}

  Users are required to specify a function of \eqn{\theta} and p and to provide starting points for the minimization,
  the design matrix (with intercept) and the response variable. Some control paramters such as, tol=1e-6,
  \eqn{\alpha=.1}, \eqn{\beta=.5}, maxit=200, maxitstart=20, cluster=NULL, display=FALSE, \eqn{\epsilon=1e-12},
  a1=.001, h1=1e-4, meth="2", lowp=.01, upp=.99, np=100 could be modified from their default.
  \eqn{\alpha} and \eqn{\beta} are parameters for line search, tol, epsilon, maxit, and a1 are parameters for quasi Newthod
  approach, maxit start is the maximum number of iteration for guessing the best start values, h1 and meth (method)
  are parameters for the gradient (method="2" is centered formula, it is possible to select "1" for right and "3"
  for five points stencil), lowp, upp and np are parameters used in the integral formula, and cluster if not NULL is a
  vector of ID to compute standard errors in longitudinal data. fun_prime_theta and fun_prime_beta are the gradient functions of the quantile function with respect to \eqn{\theta} and \eqn{\beta}
}
\value{
An object of class \dQuote{\code{niqr}}, a list containing the following items:
\item{call}{the matched call.}
\item{x}{the optimal \eqn{\theta} values.}
\item{se}{standard errors for \eqn{\theta}.}
\item{fx}{the value of the minimized integrated loss function.}
\item{gx}{the gradient calculated in the optimal points.}
\item{Hx}{the hessian of quasi newthon method}
\item{omega}{gradient matrix used in the sandwich formaula to compute standard errors.}
\item{covar}{variance covariance matrix \eqn{H^(-1)\Omega H^(-1).}}
\item{p.star.y}{the CDF calculated in the optimal x values.}
\item{code}{the convergence code, 0 for convergence, 1 for maxit reached, 2 no more step in the gradient.}
\item{internal}{a list containing some initial and control object.}
}
\references{
Frumento, P., and Bottai, M. (2015). \emph{Parametric modeling of quantile regression coefficient functions}.
Biometrics, doi: 10.1111/biom.12410.
}
\author{
Gianluca Sottile \email{gianluca.sottile@unipa.it}
}

\seealso{
\code{\link{summary.niqr}}, \code{\link{plot.niqr}}, \code{\link{predict.niqr}},
for summary, plotting, and prediction.
\code{\link{testfit.niqr}} for goodness of fit.
}
\examples{

set.seed(1234)
n <- 300
x <- runif(n)
fun <- function(theta, p){
  beta0 <- theta[1] + exp(theta[2]*p)
  beta1 <- theta[3] + theta[4]*p
  cbind(beta0, beta1)}
beta <- fun(c(1,1,1,1), runif(n))
y <- beta[, 1] + beta[, 2]*x
model <- niqr(fun=fun, x0=rep(.5, 4), X=cbind(1,x), y=y)

\dontrun{
# NOT RUN---qgamma function
set.seed(1234)
n <- 1000
x <- runif(n)
fun2 <- function(theta, p){
  beta0 <- theta[1] + qgamma(p, exp(theta[2]), exp(theta[3]))
  beta1 <- theta[4] + theta[5]*p
  cbind(beta0, beta1)
}
beta <- fun2(c(1,2,2,1,1), runif(n))
y <- beta[, 1] + beta[, 2]*x
model <- niqr(fun=fun2, x0=rep(.5, 5), X=cbind(1,x), y=y)

# NOT RUN---qbeta function
set.seed(1234)
n <- 1000
x <- runif(n)
fun3 <- function(theta, p){
  beta0 <- theta[1] + theta[2]*qbeta(p, exp(theta[3]), exp(theta[4]))
  beta1 <- theta[5] + theta[6]*p
  cbind(beta0, beta1)
}
beta <- fun3(c(1,1.5,.5,.2,1,1), runif(n))
y <- beta[, 1] + beta[, 2]*x
model <- niqr(fun=fun3, x0=rep(.5, 6), X=cbind(1,x), y=y)

# NOT RUN---qt function
set.seed(1234)
n <- 1000
x <- runif(n)
fun4 <- function(theta, p){
  beta0 <- theta[1] + exp(theta[2])*qt(p, 1+exp(theta[3]), exp(theta[4]))
  beta1 <- theta[5] + theta[6]*p
  cbind(beta0, beta1)
}
beta <- fun4(c(1,.5,.3,.2,1,1), runif(n))
y <- beta[, 1] + beta[, 2]*x
model <- niqr(fun=fun4, x0=rep(.5, 6), X=cbind(1,x), y=y)
}

# see the documentation for 'summary.piqr', and 'plot.piqr'

}
\keyword{models}
\keyword{regression}
