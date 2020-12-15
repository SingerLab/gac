#' Lookup copy number data from bins
#'
#' This function pulls the X copy number data based on coordinate information from the cnr$chromInfo
#'
#' @param cnr a cnr bundle
#'
#' @param chr a chromosome name, must match 'cnr$chromInfo$bin.chrom'
#'
#' @param start start coodinate
#'
#' @param end end coordinate
#'
#' @return
#'
#' Returns the X data for a set of chromosome coordinates
#'
#' 
#' @examples
#'
#' data(cnr)
#' coord <- data.frame(chr = 2, start = 550000, end = 600000)
#'
#' lookupCN(cnr, chr = coord$chr, start = coord$start, end = coord$end)
#'
#' @export
lookupCN <- function(cnr, chr = 12, start = 69200804, end = 69246466) {
    ## coord is a list object with elements $chr $start $end
    t(cnr$X[as.character(cnr$chromInfo$bin.chrom) == chr & cnr$chromInfo$bin.start > start & cnr$chromInfo$bin.end < end, ])
} ## lookupCN

