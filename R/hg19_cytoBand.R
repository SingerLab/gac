#' HG19 cytoband file from UCSC genome browser
#' 
#' @description
#'
#' Human reference GRCh37/hg19 cytoband file from UCSC Genome Browser
#'
#' File is used by SCclust to remove bins on the chromosome centromeres
#'
#' @format A matrix containg the cytoband coordinates
#'
#' \itemize{
#'    \item bins, in rows from j..n
#'    \item cellIDs, in columns from i..n
#' }
#' 
#' @docType data
#' @keywords datasets
#' @usage data(hg19_cytoBand)
#' @source \url{https://github.com/SingerLab/gac}
"hg19_cytoBand"
