#' Seasonal adjustment of weekly data
#'
#' Performs seasonal adjustment of weekly data. For more details on the usage of this function see the paper or the examples on Github.
#'
#' @import lubridate
#' @importFrom Hmisc yearDays
#' @importFrom stats AIC BIC lm median supsmu
#'
#' @param x Input time series as a numeric vector
#' @param dates a vector of class "Date", containing the data dates
#' @param r Defines the rate of decay of the weights. Should be between zero and one. By default is set to 0.8.
#' @param auto.ao.seacrh Boolean. Search for additive outliers
#' @param out.threshold t-stat threshold in outlier search. By default is 3.8
#' @param ao.list Vector with user specified additive outliers in a date format
#' @param my.k_l Numeric vector defining the number of yearly and monthly trigonometric variables. If NULL, is found automatically using the information criteria
#' @param H Matrix with holiday- and trading day factors
#' @param ic Information criterion used in the automatic search for the number of trigonometric regressors. There are thee options: aic, aicc and bic. By default uses aicc
#' @param method Decomposition type: additive or multiplicative
#'
#' @return sa Seasonally adjusted series
#' @return my.k_l Number of trigonometric variables used to model the seasonal pattern
#' @return sf Estimated seasonal effects
#' @return hol.factors Estimated holiday effects
#' @return out.factors Estimated outlier effects
#' @return beta Regression coefficients for the last year
#' @return m lm object. Unweighted OLS regression on the full sample
#' @author Tim Ginker
#' @export
#' @examples
#' # Not run
#' # Seasonal adjustment of weekly US gasoline production
#' \donttest{
#' data("gasoline.data")
#' res <- boiwsa(x = gasoline.data$y, dates = gasoline.data$date)
#' }
boiwsa <- function(x,
                   dates,
                   r = 0.8,
                   auto.ao.seacrh = TRUE,
                   out.threshold = 3.8,
                   ao.list = NULL,
                   my.k_l = NULL,
                   H = NULL,
                   ic = "aicc",
                   method = "additive") {
  # First run --------------------------------------

  if (method == "multiplicative") {
    x <- log(x)
  }

  # computing initial trend estimate with Friedman's SuperSmoother

  trend.init <- stats::supsmu(seq_along(x), x)$y

  y <- x - trend.init

  # looking for additive outliers

  if (auto.ao.seacrh) {
    auto.ao <- find_outliers(y = y, dates = dates, out.tolerance = out.threshold, H = H, my.AO.list = ao.list)
  } else {
    auto.ao <- NULL
  }

  if (length(auto.ao$ao) > 0) {
    ao.list <- lubridate::as_date(c(ao.list, auto.ao$ao))
  }

  # First run of the SA procedure
  # creating outlier variables
  if (is.null(ao.list)) {
    AO <- NULL

    nc.ao <- 0
  } else {
    AO <- my_ao(dates = dates, out.list = ao.list)
    nc.ao <- ncol(AO)
  }

  if (is.null(my.k_l)) {
    opt <- find_opt(y = y, dates = dates, H = H, AO = AO)

    if (ic == "aicc") {
      my.k_l <- opt$opt.aicc
    }

    if (ic == "aic") {
      my.k_l <- opt$opt.aic
    }

    if (ic == "bic") {
      my.k_l <- opt$opt.bic
    }
  }

  if (sum(my.k_l > 0)) {
    dates <- c(dates, "2024-09-18", "2024-09-25")
    X <- fourier_vars(k = my.k_l[1], l = my.k_l[2], dates = dates)
    AO <- rbind(AO, 0, 0)
    Xs <- cbind(X, H, AO)

    # Creating weights
    my.years <- unique(lubridate::year(dates))
    Wi <- array(0, dim = c(length(dates), length(dates), length(my.years)))
    w.i <- rep(0, length(my.years))

    for (i in seq_along(my.years)) {
      for (j in seq_along(my.years)) {
        w.i[j] <- r^(abs(j - i))
      }
      w <- NULL
      for (k in seq_along(my.years)) {
        w <- c(w, rep(w.i[k], sum(year(dates) == my.years[k])))
      }
      diag(Wi[, , i]) <- w / sum(w)
    }

    sf <- rep(0, length(dates))
    hol.factors <- rep(0, length(dates))
    out.factors <- rep(0, length(dates))

    for (i in seq_along(my.years)) {
      Xs_t_Wi <- t(Xs) %*% Wi[, , i]
      beta <- solve(Xs_t_Wi %*% Xs) %*% Xs_t_Wi %*% y

      n.i <- (lubridate::year(dates) == my.years[i])
      sf[n.i] <- (Xs[, 1:(length(beta) - nc.ao)] %*% beta[1:(length(beta) - nc.ao)])[n.i]

      if (!is.null(H)) {
        hol.factors[n.i] <- (Xs[, (ncol(X) + 1):(ncol(X) + ncol(H))] %*% beta[(ncol(X) + 1):(ncol(X) + ncol(H))])[n.i]
      }

      if (nc.ao > 0) {
        if (!is.null(H)) {
          out.factors[n.i] <- (Xs[, (ncol(X) + ncol(H) + 1):ncol(Xs)] %*% as.matrix(beta[(ncol(X) + ncol(H) + 1):ncol(Xs)]))[n.i]
        } else {
          out.factors[n.i] <- (as.matrix(Xs[, (ncol(X) + 1):ncol(Xs)]) %*% as.matrix(beta[(ncol(X) + 1):ncol(Xs)]))[n.i]
        }
      } else {
        out.factors <- NULL
      }
    }

    if (!is.null(out.factors)) {
      seas.out.adj <- x - sf - out.factors
    } else {
      seas.out.adj <- x - sf
    }

    # Second run --------------------------------------
    trend.init <- supsmu(seq_along(x), seas.out.adj)$y
    y <- x - trend.init

    # creating outlier variables

    if (is.null(ao.list)) {
      AO <- NULL
      nc.ao <- 0
    } else {
      AO <- my_ao(dates = dates, out.list = ao.list)
      nc.ao <- ncol(AO)
    }

    if (is.null(my.k_l)) {
      opt <- find_opt(y = y, dates = dates, H = H, AO = AO)

      if (ic == "aicc") {
        my.k_l <- opt$opt.aicc
      }

      if (ic == "aic") {
        my.k_l <- opt$opt.aic
      }

      if (ic == "bic") {
        my.k_l <- opt$opt.bic
      }
    }

    X <- fourier_vars(k = my.k_l[1], l = my.k_l[2], dates = dates)
    Xs <- cbind(X, H, AO)
    my.years <- unique(lubridate::year(dates))
    Wi <- array(0, dim = c(length(dates), length(dates), length(my.years)))
    w.i <- rep(0, length(my.years))

    for (i in seq_along(my.years)) {
      for (j in seq_along(my.years)) {
        w.i[j] <- r^(abs(j - i))
      }
      w <- NULL
      for (k in seq_along(my.years)) {
        w <- c(w, rep(w.i[k], sum(lubridate::year(dates) == my.years[k])))
      }

      diag(Wi[, , i]) <- w / sum(w)
    }

    sf <- rep(0, length(dates))
    hol.factors <- rep(0, length(dates))
    out.factors <- rep(0, length(dates))

    for (i in seq_along(my.years)) {
      Xs_t_Wi <- t(Xs) %*% Wi[, , i]
      beta <- solve(Xs_t_Wi %*% Xs) %*% Xs_t_Wi %*% y
      n.i <- (lubridate::year(dates) == my.years[i])
      sf[n.i] <- (Xs[, 1:(length(beta) - nc.ao)] %*% beta[1:(length(beta) - nc.ao)])[n.i]

      if (!is.null(H)) {
        hol.factors[n.i] <- (Xs[, (ncol(X) + 1):(ncol(X) + ncol(H))] %*% beta[(ncol(X) + 1):(ncol(X) + ncol(H))])[n.i]
      }

      if (nc.ao > 0) {
        if (!is.null(H)) {
          out.factors[n.i] <- (Xs[, (ncol(X) + ncol(H) + 1):ncol(Xs)] %*% as.matrix(beta[(ncol(X) + ncol(H) + 1):ncol(Xs)]))[n.i]
        } else {
          out.factors[n.i] <- (as.matrix(Xs[, (ncol(X) + 1):ncol(Xs)]) %*% as.matrix(beta[(ncol(X) + 1):ncol(Xs)]))[n.i]
        }
      } else {
        out.factors <- NULL
      }
    }

    if (!is.null(out.factors)) {
      seas.out.adj <- x - sf - out.factors
    } else {
      seas.out.adj <- x - sf
    }
    trend.fin <- supsmu(seq_along(x), seas.out.adj)$y

    # computing final seasonal adjusted series
    sa <- x - sf

    if (method == "multiplicative") {
      sa <- exp(sa)
      trend.fin <- exp(trend.fin)
      sf <- exp(sf)
      x <- exp(x)
    }

    lm.data <- as.data.frame(cbind(y, Xs))
    m <- lm(y ~ . - 1, data = lm.data)

    # Creating output --------------------------------------
    final_output <- list(
      sa = sa,
      my.k_l = my.k_l,
      seasonal.factors = sf,
      hol.factors = hol.factors,
      out.factors = out.factors,
      trend = trend.fin,
      beta = beta,
      m = m,
      x = x,
      dates = dates,
      ao.list = lubridate::as_date(ao.list)
    )

    class(final_output) <- "boiwsa"

    return(final_output)
  } else {
    message("Series should not be a candidate for seasonal adjustment because automatic selection found k=l=0")

    final_output <- list(
      sa = NULL,
      my.k_l = c(0, 0),
      seasonal.factors = NULL,
      hol.factors = NULL,
      out.factors = NULL,
      trend = NULL,
      beta = NULL,
      m = NULL,
      x = x,
      dates = dates,
      ao.list = NULL
    )

    class(final_output) <- "boiwsa"

    return(final_output)
  }
}
