#' Peter Andrew's arctan transformation for visualizing copy number
#'
#' The parctan transformation scales the data for visualization.  It zooms into low copy number, yet still see high valued copy number
#'
#' See http://mumdex.com
#'
#' @param x copy number data
#'
#' @examples
#' data(cnr)
#' 
#' pX <- parctan(cnr$X)
#' 
#' @export
parctan <- function(x) {

    sqrt((atan(plog(x) / plog(5))) / atan(Inf))

} ## parctan

