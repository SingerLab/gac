#' Median absolute pairwise deviance (MAPD)
#'
#' @param Xi copy number data for one cell
#'
#' @importFrom stats median
#' @importFrom stats sd
#' @importFrom utils read.delim
#' 
#' @return
#' Returns the median deviance, median deviance sd, and median deviance cv
#' 
#' @examples
#' data(cnr)
#' 
#' mapd(cnr$X[,1])
#' 
#' 
#' @export
mapd <- function(Xi) {
    az <- abs(Xi[2:length(Xi)] - Xi[1:(length(Xi)-1)])
    mz <- median(az)
    sdz <- sd(az)
    cvz <- sd(az)/mean(az)
    return(c("mapd" = mz, "mapd.sd" = sdz, "mapd.cv" = cvz))
} ## mapd

