#' Bin-to-gene index
#'
#' A table to interpolate bin copy number to individual genes.
#'
#' @description
#' 
#' This table was constructed by mapping bin regions to gene coordinates.
#'  We found the simplest way to do this, was to pull the gencode annotation
#' with biomaRt[https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html]
#' and use findOverlap with
#' GenomicRanges[https://bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html].
#' The resulting GRanges object was saved as a text file.
#' 
#'
#' @format A data frame with chromosome and position information
#' \itemize{
#'
#'   \item hgnc.symbol, HGNC gene symbol
#'   \item chrom, chromosome 
#'   \item start, gene start coordinate
#'   \item end, gene end coordinate
#'   \item width, width
#'   \item strand, gene strand
#'   \item ensembl_gene_id, ENSEMBL gene ID
#'   \item gene_biotype, gene type
#'   \item bin.id, bin ID
#' 
#' }
#' #' @docType data
#' @keywords datasets
#' @usage data(grch37.genes.5k)
#' @source \url{https://github.com/SingerLab/gac}
"grch37.genes.5k"
