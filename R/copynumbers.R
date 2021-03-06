#' Copy number quantal matrix (copynumbers)
#' 
#' @description
#' A numeric or integer matrix of bins x n.cells containing copy number quantals
#' for each bin[i] and cell[j] where bins represent a common genomic segment
#' across all cells (either a fixed with, variable binning, or .seg data).
#' A quantal is the log ratio x the ploidy.  This data can be constructed using
#' a variable length bin and CBS (Varbin algorithm) (Baslan et al 2012.), and
#' implemented on Ginkgo, or from hmmCopy.  These are upstream analyses to the
#' package.
#' 
#' The copy numbers in this alpha release is composed of uniform distributions
#' from 0 to 100 or 120 across all bins
#'
#' Alternatively, the cnq can be a binary (0,1) matrix representing presence or
#' absence of specific mutations, or a ternary (0,1,2) representing genotypes as
#' the number of alternate allele copies.
#'
#' @format A matrix containg the bin x n.cell copy number estimations
#'
#' \itemize{
#'    \item bins, in rows from j..n
#'    \item cellIDs, in columns from i..n
#' }
#'
#' 
#' @docType data
#' @keywords datasets
#' @usage data(copynumbers)
#' @source \url{https://github.com/SingerLab/GAC}
"copynumbers"
