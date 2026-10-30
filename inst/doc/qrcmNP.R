## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 4)

## ----setup--------------------------------------------------------------------
library(qrcmNP)

## ----data---------------------------------------------------------------------
set.seed(1234)
n <- 400
d <- data.frame(x1 = rexp(n), x2 = runif(n, 0, 5), x3 = rnorm(n),
                x4 = runif(n), x5 = rbinom(n, 1, 0.5))
u <- runif(n)
d$y <- 1 + qnorm(u) + (1 + 2 * u) * d$x1 + d$x2

## ----piqr---------------------------------------------------------------------
fit <- piqr(y ~ x1 + x2 + x3 + x4 + x5, formula.p = ~ slp(p, 2), data = d,
            nlambda = 40)
fit

## ----path, fig.height = 4-----------------------------------------------------
oldpar <- par(mfrow = c(1, 2))
plot(fit, xvar = "lambda")
plot(fit, xvar = "grad")
par(oldpar)

## ----gof, fig.height = 4------------------------------------------------------
best <- gof.piqr(fit, method = "BIC")
best

## ----cv, fig.height = 4-------------------------------------------------------
cv <- cv_piqr(y ~ x1 + x2 + x3 + x4 + x5, formula.p = ~ slp(p, 2), data = d,
              nfolds = 5, nlambda = 40)
cv
plot(cv)

## ----group--------------------------------------------------------------------
fitg <- piqr(y ~ x1 + x2 + x3 + x4 + x5, formula.p = ~ slp(p, 2), data = d,
             nlambda = 40, penalty = "group")
bestg <- gof.piqr(fitg, plot = FALSE)
bestg

## ----refit--------------------------------------------------------------------
m <- refit_piqr(fit, pos.lambda = best$posMinLambda)
summary(m)
summary(m, p = c(0.1, 0.5, 0.9))

## ----refit-plot, fig.height = 3.5---------------------------------------------
oldpar <- par(mfrow = c(1, 3))
plot(m, which = 1:3, ask = FALSE)
par(oldpar)

## ----censored-----------------------------------------------------------------
cens <- rexp(n, 0.05)
d$time <- pmin(d$y + 5, cens)
d$status <- as.numeric(d$y + 5 <= cens)
fitc <- piqr(survival::Surv(time, status) ~ x1 + x2 + x3 + x4 + x5,
             formula.p = ~ slp(p, 2), data = d, nlambda = 20)
gof.piqr(fitc, plot = FALSE)

## ----niqr1--------------------------------------------------------------------
set.seed(1)
n <- 400
x <- runif(n)
z <- qnorm(runif(n))
y <- 2 + 0.5 * z + 0.1 * z^3 + (1 + pnorm(z)) * x
dn <- data.frame(y, x)

# Model 1: normal intercept
m1 <- niqr(y ~ x, data = dn, x0 = c(2, 0, 1, 0),
           fun = function(theta, p) cbind(theta[1] + exp(theta[2]) * qnorm(p),
                                          theta[3] + theta[4] * p))
# Model 2: a cubic term accounts for the kurtosis (exp() keeps it positive,
# so that the quantile function is increasing in p)
m2 <- niqr(y ~ x, data = dn, x0 = c(2, 0, -2, 1, 0),
           fun = function(theta, p) {
             z <- qnorm(p)
             cbind(theta[1] + exp(theta[2]) * z + exp(theta[3]) * z^3, theta[4] + theta[5] * p)
           })
gof.niqr(m1, m2)
summary(m2)

## ----niqr-plot, fig.height = 3.5----------------------------------------------
oldpar <- par(mfrow = c(1, 3))
plot(m2, which = 1:3, ask = FALSE)
par(oldpar)

## ----testfit------------------------------------------------------------------
test.fit(m1, R = 30)
test.fit(m2, R = 30)

## ----hill---------------------------------------------------------------------
set.seed(2)
dose <- rep(10^seq(-2, 2, length.out = 20), each = 15)
fun.h <- function(theta, p) cbind(theta[1] + exp(theta[2]) * qnorm(p), theta[3], theta[4])
fun2.h <- function(beta, X) beta[, 1] + (beta[, 2] - beta[, 1]) / (1 + 10^(beta[, 3] - log10(X[, 2])))
yh <- fun2.h(fun.h(c(1, log(0.5), 10, 0), runif(length(dose))), cbind(1, dose))
mh <- niqr(yh ~ dose, fun = fun.h, fun2 = fun2.h, x0 = c(0, 0, 5, 0.5))
summary(mh)

## ----hill-plot, fig.height = 4------------------------------------------------
dd <- data.frame(dose = 10^seq(-2, 2, length.out = 100))
qf <- predict(mh, type = "QF", p = c(0.05, 0.5, 0.95), newdata = dd, se = FALSE)
plot(dose, yh, log = "x", pch = 20, col = "gray", xlab = "dose", ylab = "response")
matlines(dd$dose, qf, lty = c(2, 1, 2), col = 1)

