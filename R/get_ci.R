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
  stopifnot(which %in% c("high", "low"))
  if (isTRUE(parametric)) {
    switch(which,
      high = mean(data) + sd(data) * sd_mult,
      low = mean(data) - sd(data) * sd_mult
    )
  } else {
    switch(which,
      high = quantile(data, probs = (1 + ci) / 2, type = 1),
      low = quantile(data, probs = (1 - ci) / 2, type = 1)
    )
  }
}
