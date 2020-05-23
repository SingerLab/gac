#' Chromosome coordinates and additional data of the bins
#' 
#' @description
#' A data.frame containing the chromsome and end coordiantes for all the bins
#' of the single-cell copy number matrix `copynumbers`.  This table is used for
#' SCclust, some data may be recurrent now to preserve compatibility with that
#' package and the functions here.
#' This table is cross referenced to create the `genes` table, an interpolated
#' table to link the bins to genes
#' 
#' @format A data frame with chromosome and position information
#' \itemize{
#' 
#'  \item bin.chrom, object class `factor` corresponding to the chromosome
#' of the bin.
#'  \item bin.start, starting coordinate of the bin
#'  \item bin.start.abspos, cumulative start positions of the bins
#'  \item bin.end, end coordinate of the bin
#'  \item bin.length, size of the bin
#'  \item mappable.positions, number of mappable positions
#'  \item gc.content, percent gc content
#'  \item chrom.numeric, a numeric representation of the autosomes (e.g. 1..24 for hg19, 1..23 for mm10, 1..30 for bta9
#'  \item chrom, alternate representation of the chromosome to bind to X for SCclust
#'  \item chrompos, alternate representation of the bin.start to bind to X
#'   for SCclust
#'  \item abspos, alternate representation of the chromosome to bind to X for
#'   SCclust
#' }
#' 
#' @docType data
#' @keywords datasets
#' @usage data(chromInfo)
#' @source \url{https://github.com/SingerLab/gac}
"chromInfo"
