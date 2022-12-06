#' Correlation to dataframe
#' @description Helper function for bootstrap_traits_multivariate that
#'  extracts results from a correlation matrix
#' @param corr correlation matrix
#' @return A data.frame of correlations
#' @references Modified from https://stackoverflow.com/a/23476844/2055765
#' @export
#' @importFrom stats as.dist
#' @examples
#' x <- matrix(ncol = 5, rnorm(20))
#' colnames(x) <- letters[1:5]
#' cor(x) |> cor_to_df()
#'
cor_to_df <- function(corr) {
  if (nrow(corr) != ncol(corr)) stop("needs square matrix")
  corr <- as.dist(corr)
  A <- attr(corr, "Size")
  B <- attr(corr, "Labels")
  data.frame(
    row = B[unlist(lapply(2:A, function(x) x:A))],
    col = rep(B[-length(B)], (length(B) - 1):1),
    value = as.vector(corr)
  )
}
