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
  size <- attr(corr, "Size")
  labels <- attr(corr, "Labels")
  data.frame(
    row = labels[unlist(lapply(2:size, function(x) x:size))],
    col = rep(labels[-length(labels)], (length(labels) - 1):1),
    value = as.vector(corr)
  )
}
