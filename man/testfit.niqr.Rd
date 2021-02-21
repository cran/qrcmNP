\name{testfit.niqr}
\alias{testfit.niqr}
\title{
Goodness-of-Fit Test
}
\description{
Goodness-of-fit test for a model
fitted with \code{\link{niqr}}. The Kolmogorov-Smirnov statistic and the Cramer-Von Mises statistic
are computed. Their distribution under the null hypothesis is estimated
with Monte Carlo (see \sQuote{Details}).
}
\usage{
testfit.niqr(obj, R = 100)
}
\arguments{
  \item{obj}{
    an object of class \dQuote{\code{niqr}}.
  }
  \item{R}{
   number of Monte Carlo replications.
}
}
\details{
  This function permits assessing goodness of fit by testing the null hypothesis
that the \acronym{CDF} values follow a \eqn{U(0,1)} distribution, indicating that
the model is correctly specified.
Since the \acronym{CDF} values depend on estimated parameters, the distribution of
the test statistic is not known. To evaluate it, the model is fitted on \kbd{R} simulated datasets
generated under the null hypothesis.
}
\value{
a matrix with columns \code{statistic} and \code{p.value},
reporting the Kolmogorov-Smirnov and Cramer-Von Mises statistic and the associated
p-values evaluated with Monte Carlo.
}
\references{
Frumento, P., and Bottai, M. (2015). \emph{Parametric modeling of quantile regression coefficient functions}.
Biometrics, doi: 10.1111/biom.12410.
}
\author{
Gianluca Sottile \email{gianluca.sottile@unipa.it}
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
\dontrun{testfit.niqr(model, R=100)}
}
\keyword{htest}
