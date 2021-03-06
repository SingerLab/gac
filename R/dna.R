#'  dna: a cnr with bulk DNA like copy number data
#'
#' @description
#' 
#' The cnr object is a list of four relational matrices. The bins, genes,
#' annotation, qc, chromInfo, and gene.index.  The structure is inspired by
#' Scanpy's AnnData to which cleverly integrates complex data into a simple
#' architecture.
#'
#'
#' @format An object class list containg a rounded CNR
#'
#' \itemize{
#' 
#'   \item X, An integer matrix of bins x n.cells containing DNA copy number ratio
#' estimations for each bin[i] and cell[j] Where bins represent a common genomic
#' segment across all cells (either a fixed with, variable binning, or .seg data).  This data can be constructed using a variable length bin and CBS (Varbin
#' algorithm) (Baslan et al 2012.), and implemented on Ginkgo, or from hmmCopy.
#' These are upstream analyses to the package.
#'
#' For mutations in single-cells, X can be a binary incidence (0,1) matrix
#' representing presence or absence of specific mutations, or a ternary (0,1,2)
#' representing genotypes as the number of alternate allele copies
#'
#'   \item genes, gene copy number interpolation from bins.  The genes matrix
#' is an interpolated, transposed, expansion of bins. The expansion is
#' constructed internally using the expand2genes function.
#'
#'   \item Y, phenotypic data of single-cells, contains cells as rows, and
#' phentypes in columns. Phenotypes can be information about individual samples,
#' or if same-cell methods were used, the RNA expression from the same cell.
#'#' 
#'   \item qc, quality control metrics. This matrix contains additional metadata
#' that is technical, e.g. number of reads, MAPD estimates, and the PASS/FAIL
#' qc.status for individual cells.  contains cells as rows and metadata as columns
#'
#'  \item chromInfo, ordered chromosome information for the bins
#'
#'  \item gene.index, table to map bins to genes
#' #'   ...
#' }
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @usage data(dna)
#' 
#' @source \url{https://github.com/SingerLab/gac"}
"dna"
