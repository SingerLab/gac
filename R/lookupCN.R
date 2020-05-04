#' Lookup copy number data from bins
#'
#' This function pulls the X copy number data based on coordinate information from the cnr$chromInfo
#'
#' @param cnr a cnr bundle
#'
#' @param coord a `list` with elements named chr start end
#'
#' @return
#'
#' Returns the X data for a set of chromosome coordinates
#'
#' 
#' @examples
#'
#' data(cnr)
#' coords <- list(chr = 2, start = 550000, end = 600000)
#'
#' lookupCN(cnr, coord = coords)
#' 
#' @export
lookupCN <- function(cnr, coord) {
    ## coord is a list object with elements $chr $start $end
    t(cnr$X[cnr$chromInfo$chr == coord$chr & cnr$chromInfo$end > coord$start & cnr$chromInfo$end < coord$end, ])
} ## lookupCN

