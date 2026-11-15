#' Get different distributions
#' @description Internal helper function for drawing traits
#' @param parm1 the first parameter in a fitted distribution
#' @param parm2 the second parameter in a fitted distribution
#' @param n Number of draws.  Defaults to 1
#' @param type Distribution type to fit.
#' One of "normal", "lognormal", or "beta"
#' @returns numeric vector of length `n` from required distribution.
#' @importFrom stats rnorm rlnorm rbeta
#' @importFrom glue glue
#' @keywords internal

distribution_handler <- function(parm1, parm2, n = 1, type) {
  switch(type,
    normal = rnorm(n = n, mean = parm1, sd = parm2),
    lognormal = rlnorm(n = n, meanlog = parm1, sdlog = parm2),
    beta = rbeta(n = n, shape1 = parm1, shape2 = parm2),
    stop(glue("Type must be 'normal', 'lognormal', or 'beta', not {{type}}"))
  )
}
