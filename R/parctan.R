#' Peter Andrew's arctan transformation for visualizing copy number
#'
#' The parctan transformation scales the data for visualization.  It zooms into low copy 
#' number, yet still see high valued copy number
#'
#' See http://mumdex.com
#'
#' @param x copy number data
#'
#' @return
#'
#' Returns a copy number object on Peter Andrew's arctan transformation
#' 
#' @examples
#' data(cnr)
#' 
#' pX <- parctan(cnr$X)
#' 
#' @export
parctan <- function(x) {

    sqrt((atan(plog(x) / plog(5))) / atan(Inf))

} ## end parctan

#' internal
#'
#' calculates a plog for parctan()
#'
#' @param x vector or matrix of copy numbers
#'
#' @return
#' Returns plog transformed data
#' 
#' @keywords internal
plog <- function(x) {
    log10(1 + x^2)
}

