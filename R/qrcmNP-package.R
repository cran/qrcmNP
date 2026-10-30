#' qrcmNP: Nonlinear and Penalized Quantile Regression Coefficients Modeling
#'
#' Quantile regression coefficients modeling (QRCM, Frumento and Bottai, 2016)
#' describes the regression coefficients \eqn{\beta(p)} of a quantile
#' regression model as parametric functions of the quantile order \eqn{p}, and
#' estimates all quantiles simultaneously by minimizing an integrated loss
#' function. 'qrcmNP' extends this framework in two directions:
#'
#' * [piqr()]: penalized QRCM for covariate selection (Sottile et al., 2020),
#'   with lasso and group lasso penalties, tuning by AIC/BIC ([gof.piqr()]) or
#'   cross-validation ([cv_piqr()]), and post-selection inference
#'   ([refit_piqr()]);
#' * [niqr()]: nonlinear parametric quantile models, in which the quantile
#'   function can be a nonlinear function of the parameters, of the covariates,
#'   or both (Bottai and Cilluffo, 2020), with sandwich, cluster-robust and
#'   bootstrap standard errors, information criteria ([gof.niqr()]) and a
#'   goodness-of-fit test ([test.fit.niqr()]).
#'
#' @references
#' Frumento, P. and Bottai, M. (2016). Parametric modeling of quantile
#' regression coefficient functions. *Biometrics*, 72(1), 74-84.
#' \doi{10.1111/biom.12410}
#'
#' Sottile, G., Frumento, P., Chiodi, M. and Bottai, M. (2020). A penalized
#' approach to covariate selection through quantile regression coefficient
#' models. *Statistical Modelling*, 20(4), 369-385.
#' \doi{10.1177/1471082X19825523}
#'
#' Bottai, M. and Cilluffo, G. (2020). Nonlinear parametric quantile models.
#' *Statistical Methods in Medical Research*, 29(12), 3757-3769.
#' \doi{10.1177/0962280220941159}
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom graphics abline axis lines matplot par plot points polygon text
#' @importFrom grDevices adjustcolor gray.colors n2mfrow
#' @importFrom stats approx coef delete.response model.frame model.matrix
#'   model.response model.weights nobs pnorm ppoints predict printCoefmat
#'   runif sd terms vcov .getXlevels
#' @importFrom utils menu
#' @importFrom qrcm test.fit
NULL

#' @importFrom qrcm slp
#' @export
qrcm::slp

#' @importFrom qrcm plf
#' @export
qrcm::plf

#' @export
qrcm::test.fit
