#' internal
#'
#' calculates a plog for parctan()
#'
#' @param x vector or matrix of copy numbers
#'
#' @return
#' Returns plog transformed data
#'
#' @examples
#' data(cnr)
#' 
#' plog(cnr$X[,1])
#' 
#' @export
plog <- function(x) {
    log10(1 + x^2)
}

