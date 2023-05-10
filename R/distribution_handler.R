#' Get different distributions
#' @description Internal helper function for drawing traits
#' @param parm1 the first parameter in a fitted distribution
#' @param parm2 the second parameter in a fitted distribution
#' @param n Number of draws.  Defaults to 1
#' @param type Distribution type to fit.
#' One of "normal", "lognormal", or "beta"
#' @returns numeric vector of length `n` from required distribution.
#' @importFrom stats rnorm rlnorm rbeta
#' @keywords internal

distribution_handler <- function(parm1, parm2, n = 1, type) {
  if (type == "normal") {
    return(rnorm(n = n, mean = parm1, sd = parm2))
  }
  if (type == "lognormal") {
    return(rlnorm(n = n, meanlog = parm1, sdlog = parm2))
  }
  if (type == "beta") {
    return(rbeta(n = n, shape1 = parm1, shape2 = parm2))
  }
}
