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


  if (distribution_type == "normal") {
    try(fit <- fitdist(data = data, distr = "norm",
                        keepdata = F, method = "mme"))
    if (!exists("fit")) {
      fit <- fitdist(data = data, distr = "norm", keepdata = F, method = "mle")
      }

    gof <- gofstat(fit)
    return(data.frame(distribution_type = "normal",
                      parm1 = fit$estimate[[1]],
                      parm2 = fit$estimate[[2]],
                      ks = gof$ks,
                      cvm = gof$cvm,
                      ad = gof$ad)
    )


  }


  if (distribution_type == "lognormal") {
    try(fit <- fitdist(data = data,
                        distr = "lnorm", keepdata = F, method = "mme"))
    if (!exists("fit")) {
      fit <- fitdist(data = data,
                     distr = "lnorm", keepdata = F, method = "mle")
      }

    gof <- gofstat(fit)
    return(data.frame(distribution_type = "lognormal",
                      parm1 = fit$estimate[[1]],
                      parm2 = fit$estimate[[2]],
                      ks = gof$ks,
                      cvm = gof$cvm,
                      ad = gof$ad)
    )


  }


  if (distribution_type == "beta") {
    try(fit <- fitdist(data = data,
                       distr = "beta", keepdata = F, method = "mme"))
    if (!exists("fit")) {
      fit <- fitdist(data = data,
                     distr = "beta", keepdata = F, method = "mle")
    }
    gof <- gofstat(fit)
    return(data.frame(distribution_type = "beta",
                      parm1 = fit$estimate[[1]],
                      parm2 = fit$estimate[[2]],
                      ks = gof$ks,
                      cvm = gof$cvm,
                      ad = gof$ad)
    )


  }


}
