#' Do something quietly
#' @description Internal helper function
#' used to prevent printed messages
#' @param x Code to be executed quietly
#' @keywords internal

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  force(x)
}
