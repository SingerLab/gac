#' Bin chromosome-basepair coordinates
#' 
#' @description
#' A dataset containing the chromsome and end coordiantes for all the bins of the single-cell copy number matrix `copynumbers`
#' This table is cross referenced to create the `genes` table, an interpolated table to link the bins to genes
#'
#' @format A data frame with chromosome and position information
#' \itemize{
#'
#'  \item chr, chromosome
#' 
#'  \item end, base pair coordinate
#' 
#' }
#' @docType data
#' @keywords datasets
#' @usage data(chromInfo)
#' @source \url{https://github.com/SingerLab/GAC}
"chromInfo"

