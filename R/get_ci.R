#' Get confidence intervals
#' @description Internal helper function for generating CIs.
#' @param data Numeric vector
#' @param sd_mult SD multiplier.  defaults to 1
#' @param ci 0 to 1
#' @param which "high" or "low"
#' @param parametric if TRUE, sd used, if FALSE, quantile
#' @returns numeric vector length 1
#' @importFrom stats quantile
#' @keywords internal

get_ci <- function(data, sd_mult = 1, ci = 0.95, which, parametric = TRUE) {
  if (isTRUE(parametric)) {
    if (which == "high") {
      return(mean(data) + sd(data) * sd_mult)
    }
    if (which == "low") {
      return(mean(data) - sd(data) * sd_mult)
    }
  } else {
    if (which == "high") {
      return(quantile(data, probs = (1 + ci) / 2, type = 1))
    }
    if (which == "low") {
      return(quantile(data, probs = (1 - ci) / 2, type = 1))
    }
  }
}
