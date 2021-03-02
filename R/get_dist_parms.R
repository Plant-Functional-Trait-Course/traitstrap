#' Get distribution parameters
#' @description Internal helper function for
#' calculating parametric distribution parameters.
#' @param data Numeric vector
#' @param distribution_type Character.
#' One of either "normal", "lognormal", or "beta"
#' @note This function will first try fitting using mme, then mle
#' @importFrom fitdistrplus fitdist gofstat
#' @keywords internal

get_dist_parms <- function(data, distribution_type) {

  if (length(distribution_type) > 1) {
    stop("This function only handles one distribution type at a time.")

  }

  # Check that distribution type is supported
  if (!distribution_type %in% c("normal", "beta", "lognormal")) {
    stop("Unsupported distribution type")
  }

  #Set distr for fitdist
  if (distribution_type == "normal") {
    distr <- "norm"
  }

  if (distribution_type == c("lognormal")) {
    distr <- "lnorm"
  }

  if (distribution_type == "beta") {
    distr <- "beta"
  }

  #Hack to have fit a distribution with only one point for normal or lognormal
  if (distribution_type %in% c("normal", "lognormal") & length(data) == 1) {
    data <- replicate(n = 2, expr = data)

  }

  #Fit model
  fit <- tryCatch({
    fitdist(data = data,
            distr = distr,
            keepdata = FALSE, method = "mme")
  }, error = function(e) {
    return(fitdist(data = data, distr = distr,
                          keepdata = FALSE, method = "mle"))

    }
           )

  #Goodness of fit test
  gof <- gofstat(fit)

  #Ensure NULL values for sd don't break anything
  if (is.null(fit$sd)) {
    fit$sd <- c(NA, NA)
  }

  #Return results
  return(data.frame(distribution_type = distribution_type,
                    parm1 = fit$estimate[[1]],
                    parm2 = fit$estimate[[2]],
                    sd1 = fit$sd[[1]],
                    sd2 = fit$sd[[2]],
                    ks = gof$ks,
                    cvm = gof$cvm,
                    ad = gof$ad)
         )


}
