#' Find optimal number of fourier variables
#'
#' Searches through the model space to identify the best number of trigonometric variables, with the lowest AIC, AICc or BIC value.
#'
#' @importFrom stats supsmu lm AIC BIC
#' @import lubridate
#'
#' @param x Numeric vector. Time series to seasonally adjust
#' @param dates a vector of class "Date", containing the data dates
#' @param H (optional) Matrix with holiday and trading day variables
#' @param AO (optional) Matrix with additive outlier variables
#' @param method Decomposition method: "additive" or "multiplicative". By default uses the additive method
#' @param l.max maximal number of the monthly cycle variables to search for
#' @param k.max maximal number of the yearly cycle variables to search for
#' @param by step size in the search
#'
#' @return list with the optimal number of (yearly and monthly) fourier variables according to AIC, AICc and BIC
#' @export
#'
#' @examples
#'
#' data(gasoline.data)
#'
#' res <- find_opt(x = gasoline.data$y, dates = gasoline.data$date)
#' print(res)
#'
find_opt <- function(y, dates, H = NULL, AO = NULL) {
  # y - detrended dependent variable
  # H - holiday and trading day effects (as matrix)
  # AO - additive outlier variables (as matrix)

  aic0 <- matrix(NA, nrow = length(seq(6, 42, 6)), ncol = length(seq(6, 42, 6)))
  aicc0 <- matrix(NA, nrow = length(seq(6, 42, 6)), ncol = length(seq(6, 42, 6)))
  bic0 <- matrix(NA, nrow = length(seq(6, 42, 6)), ncol = length(seq(6, 42, 6)))

  for (i in seq_along(seq(6, 42, 6))) {
    for (j in seq_along(seq(6, 42, 6))) {
      X <- fourier_vars(k = (i - 1) * 6, l = (j - 1) * 6, dates)

      X <- cbind(X, H, AO)

      if (is.null(X)) {
        m <- stats::lm(y ~ -1)
      } else {
        m <- stats::lm(y ~ X - 1)
      }

      aic0[i, j] <- stats::AIC(m)
      aicc0[i, j] <- stats::AIC(m) + 2 * length(m$coefficients) * (length(m$coefficients) + 1) / (length(m$residuals) - length(m$coefficients) + 1)
      bic0[i, j] <- stats::BIC(m)
    }
  }

  opt.aic <- (which(aic0 == min(aic0), arr.ind = TRUE) - 1) * 6 # optimal number of terms
  opt.aicc <- (which(aicc0 == min(aicc0), arr.ind = TRUE) - 1) * 6
  opt.bic <- (which(bic0 == min(bic0), arr.ind = TRUE) - 1) * 6

  return(list(opt.aic = opt.aic, opt.aicc = opt.aicc, opt.bic = opt.bic))
}
