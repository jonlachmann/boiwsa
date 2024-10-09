#' Find additive outliers
#'
#' Searches for additive outliers using the method described in Appendix C of Findley et al. (1998).
#' If the number of trigonometric variables is not specified will search automatically through the model space to identify the best number of trigonometric variables, with the lowest AIC, AICc or BIC value.
#'
#' @importFrom stats AIC BIC lm median supsmu
#' @import lubridate
#'
#' @param x Numeric vector. Time series to seasonally adjust
#' @param dates a vector of class "Date", containing the data dates
#' @param out.tolerance t-stat threshold for outliers (see Findley et al., 1998)
#' @param my.AO.list (optional) Vector with user defined additive outlier variables
#' @param H (optional) Matrix with holiday and trading day variables
#' @param my.k_l (optional) Vector with the number of fourier terms to capture the yearly and monthly cycle. If NULL, would perform automatic search using AICc criterion
#' @param method Decomposition method: "additive" or "multiplicative". By default uses the additive method
#'
#' @return my.k_l
#' @return ao list of AO dates
#' @references Findley, D.F., Monsell, B.C., Bell, W.R., Otto, M.C. and B.C Chen (1998). New capabilities and methods of the X-12-ARIMA seasonal-adjustment program. Journal of Business & Economic Statistics, 16(2), pp.127-152.
#' @export
#'
#' @examples
#' \donttest{
#' # Not run:
#' # Searching for additive outliers in Gasoline data
#' data(gasoline.data)
#' ao_list <- find_outliers(x = gasoline.data$y, dates = gasoline.data$date)
#' }
#'
find_outliers <- function(y,
                          dates,
                          out.tolerance = out.threshold,
                          my.AO.list = NULL,
                          H = NULL,
                          my.k_l = NULL) {
  # y -detrended variable
  # out.tolerance - t-stat threshold
  # predefined additive outlier variables
  # my.k_l - number of yearly and monthly fourier variables

  if (is.null(my.k_l)) {
    if (is.null(my.AO.list)) {
      AO <- NULL
    } else {
      AO <- my_ao(dates = dates, out.list = my.AO.list)
    }
    opt <- find_opt(y = y, dates = dates, H = H, AO = AO)
    my.k_l <- opt$opt.aicc
  }

  if (sum(opt$opt.aicc) > 0) {
    X <- fourier_vars(k = my.k_l[1], l = my.k_l[2], dates = dates)
    Xs <- cbind(X, H, AO)
    err <- y - Xs %*% solve(t(Xs) %*% Xs) %*% t(Xs) %*% y
    sig_R <- 1.49 * stats::median(abs(err))
    f.sel.pos <- NULL
    out.search.points <- as.integer((seq_along(dates))[!dates %in% my.AO.list])
    run <- TRUE

    Xs_t <- t(Xs)
    while (run) {
      cat(ncol(Xs), "")
      Ts <- numeric(length(out.search.points))
      ts_idx <- 1
      Xst2_inv <- solve(crossprod(Xs))
      Xst_y <- t(Xs) %*% y
      for (t in out.search.points) {
        AOt <- rep(0, length(dates))
        AOt[t] <- 1
        Xst2_inv_t <- rankUpdateInverse(Xst2_inv, Xs_t, AOt)
        Xst_y_t <- rbind(Xst_y, t(AOt) %*% y)
        Tt <- (Xst2_inv_t %*% Xst_y_t)[ncol(Xs) + 1] / (diag(Xst2_inv_t * sig_R^2)[ncol(Xs) + 1]^0.5)
        Ts[ts_idx] <- abs(Tt)
        ts_idx <- ts_idx + 1
      }
      if (max(Ts) >= out.tolerance) {
        AOt <- rep(0, length(dates))
        AOt[out.search.points[which.max(Ts)]] <- 1

        f.sel.pos <- c(f.sel.pos, out.search.points[which.max(Ts)])
        out.search.points <- out.search.points[-which.max(Ts)]

        Xs <- cbind(Xs, AOt)
        Xs_t <- t(Xs)
      }

      if (max(Ts) < out.tolerance) {
        run <- FALSE
      }
    }

    # Backward deletion
    if (length(f.sel.pos) > 0) {
      run <- TRUE
      f.sel.ao.dates <- dates[f.sel.pos]
    } else {
      f.sel.ao.dates <- NULL
    }

    cat("\n")
    while (run) {
      cat(".")
      AObd <- my_ao(dates = dates, out.list = lubridate::as_date(c(my.AO.list, f.sel.ao.dates)))

      Xst <- cbind(X, H, AObd)
      Xst_t <- t(Xst)
      Xst2_inv <- solve(Xst_t %*% Xst)

      err <- y - Xst %*% Xst2_inv %*% Xst_t %*% y
      sig_R <- 1.49 * stats::median(abs(err))
      Tt <- abs((Xst2_inv %*% Xst_t %*% y) / (diag(Xst2_inv * sig_R^2)^0.5))[(ncol(Xst) - length(f.sel.ao.dates) + 1):ncol(Xst)]

      if (min(Tt) < out.tolerance) {
        f.sel.ao.dates <- f.sel.ao.dates[-which.min(Tt)]
      } else {
        run <- FALSE
      }

      if (length(f.sel.ao.dates) == 0) {
        run <- FALSE
      }
    }

    if (length(f.sel.ao.dates) == 0) {
      f.sel.ao.dates <- NULL
    } else {
      f.sel.ao.dates <- f.sel.ao.dates[order(f.sel.ao.dates)]
    }

    return(list(ao = f.sel.ao.dates, my.k_l = my.k_l))
  } else {
    return(list(ao = NULL, my.k_l = my.k_l))
  }
}

#' Update the inverse of a cross product of a matrix X when adding a new column v
#' @param X_inv The inverse (X^T X)^-1 before adding the new column
#' @param X_t The transpose of X, i.e. X^T
#' @param v The column to add
#' @return The inverse of ([X v]^T [X v])^-1
rankUpdateInverse <- function(X_inv, X_t, v) {
  u1 <- X_t %*% v
  u2 <- X_inv %*% u1
  d <- as.numeric(1 / (t(v) %*% v - t(u1) %*% u2))
  u3 <- d * u2
  F11_inv <- X_inv + d * u2 %*% t(u2)
  XtX_inv <- rbind(cbind(F11_inv, -u3), c(-u3, d))
  return(XtX_inv)
}
