#' Single-cell phenotype annotation
#'
#' A dataset containing single-cell phenotypic data
#' 
#'
#' @format A data frame or matrix n.cell rows x y.phenotypes columns
#'
#' The phenotype table contains several traits that mimic what you would collect as clinical phenotypes.  Binary e.g. yes/no codes a 0/1. Categorical, usually strings, e.g. low, medium, high.  Quantitative numerical measurments from assays, or as simple as height or body weight.  A few random traits to use as controls.
#'
#' \itemize{
#'    \item cellID, single-cell ID. By default, the rownames will be used by
#' buildCNR as the cellID. It is imperative that the rownames are used as cellID,
#' if you have a column called `cellID` it will be re-written. 
#'
#'   \item binary1, binary trait 1
#'   \item binary2, binary trait 2
#'
#'   \item category1, categorical trait 1
#'   \item category2, categorical trait 2
#'
#'   \item quantitative1, quantitative trait 1
#'   \item quantitative2, quantitative trait 2
#'
#'   \item random1, random trait 1
#'   \item random2, random trait 2
#' }
#' 
#' @docType data
#' @keywords datasets
#' @usage data(pheno)
#' @source \url{https://github.com/SingerLab/gac/}
"pheno"
