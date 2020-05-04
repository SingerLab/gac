#' Lookup copy number data from bins
#'
#' This function pulls the set of copy number data for a specific set of coordinates from a table, e.g. the CNR$chromInfo, 
#'
#' @param df `data.frame` containing chrom, and chrompos
#'
#' @param coord a `list` with elements named chr start end
#'
#' @examples
#'
#'
#' data(cnr)
#' coords <- list(chr = 2, start = 550000, end = 600000)
#'
#' lookupCN(df = cnr$chromInfo, coord = coords)
#' 
#' @export
lookupCN <- function(df, coord) {
    ## coord is a list object with elements $chr $start $end
    t(df[df$chrom == coord$chr & df$chrompos > coord$start & df$chrompos < coord$end, 4:ncol(df)])
} ## lookupCN

