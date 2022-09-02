#' Subset a set of bins or genes from a CNR object
#'
#' This function pulls out a section of the genome for further analysis. 
#'
#' @param cnr a cnr bundle
#'
#' @param bins a vector of bins to keep
#'
#' @param genes a vector of genes to keep
#'
#' @param chrom a vector of chromosomes, or single chromosome if genomic coordinates are used
#'
#' @param start a start position
#'
#' @param end an end position
#'
#' @param ... additional parameters e.g. padding to add bins when coordinates are used
#' 
#' @return
#'
#' Returns a cnr object with only the desired bins.  This function also subsets the gene index and chromInfo as well
#'
#' @examples
#' data(cnr)
#' data(segCol)
#'
#' ## subset based on bins
#' bins1.100 <- subsetCNR(cnr, bins = 1:100)
#'
#' ## subset based on genes
#' genes.3 <- subsetCNR(cnr, genes = c("CDK4", "MDM2", "HMGA2", "TP53"))
#'
#' ## subset based on chromosomes
#' chrom6.12 <- subsetCNR(cnr, chrom = c(6, 12))
#' 
#' ## subset based on genomic coordinates
#' chrom12q13.15 <- subsetCNR(cnr, chrom = 12, start = 46000000, end = 71000000)
#'
#' @importFrom assertthat assert_that
#' 
#' @export
subsetCNR <- function(cnr, bins = NULL, genes = NULL, chrom = NULL, start = NULL, end = NULL, ...) {

    if(!is.null(bins)) {
        assertthat::assert_that(is.null(genes))
        assertthat::assert_that(is.null(chrom))
        assertthat::assert_that(is.null(start))
        assertthat::assert_that(is.null(end))
        
        cnr <- subset_on_bins(cnr, bins = bins)
    }
    
    if(!is.null(genes)) {
        assertthat::assert_that(is.null(bins))
        assertthat::assert_that(is.null(chrom))
        assertthat::assert_that(is.null(start))
        assertthat::assert_that(is.null(end))
        
        cnr <- subset_on_genes(cnr, genes = genes, ...)
    }

    if(!is.null(chrom) & is.null(start) & is.null(end)) {
        assertthat::assert_that(is.null(bins))
        assertthat::assert_that(is.null(genes))

        cnr <- subset_on_chrom(cnr, chrom = chrom)
    }

    if(!is.null(chrom) & !is.null(start) & !is.null(end)) {
        assertthat::assert_that(is.null(bins))
        assertthat::assert_that(is.null(genes))

        cnr <- subset_on_coordinates(cnr, chrom = chrom,
                                     start = start, end = end, ...)
    }
    
    return(cnr)
    
}


#' subset based on bins
#' @param cnr a cnr bundle
#'
#' @param bins a list of bins to subset
#' 
#' @keywords internal
#' @noRd
subset_on_bins <- function(cnr, bins) {

    ## subset X based on bins
    nX <- cnr$X[bins, ]

    ## subset gene index based on bins
    g2 <- cnr$gene.index[cnr$gene.index$bin.id %in% bins, ]
    rownames(g2) <- g2$hgnc.symbol
    
    ## subset genes matrix based on bins
    gg <- gsub("-", ".", g2$hgnc.symbol)
    nGenes <- cnr$genes[, gg]
        
    ## subset chromInfo based on bins
    nCI <- cnr$chromInfo[bins,]

    ## if expression is present
    if(!is.null(cnr[["expr"]])) {
        nE <- cnr$expr[, gg]
    }

    ## if DDRC.df is present
    if(!is.null(cnr[["DDRC.df"]])) {
        nDDRC.df <- cnr$DDRC.df[bins, ]
        nDDRC.g <- cnr$DDRC.g[, gg]
    }
    
    cnr[["X"]] <- nX
    cnr[["genes"]] <- nGenes
    cnr[["gene.index"]] <- g2
    cnr[["chromInfo"]] <- nCI
    if(!is.null(cnr[["expr"]])) {
        cnr[["expr"]] <- nE
    }
   
    if(!is.null(cnr[["DDRC.df"]])) {
       cnr[["DDRC.df"]] <- nDDRC.df
       cnr[["DDRC.g"]] <- nDDRC.g
    }
    
    return(cnr)
} # end subset_on_bins


#' subset based on genes
#' @param cnr a cnr bundle
#'
#' @param genes a list of genes to subset
#'
#' @param all return all genes within on the bins of the genes of interst
#' 
#' @keywords internal
#' @noRd
subset_on_genes <- function(cnr, genes, all = TRUE) {

    ## transforming genes to column names
    gg <- gsub("-", ".", genes)
    
    ## subset genes matrix based on genes
    b2 <- as.numeric(gac::mark.genes(cnr, gene.list = genes))

    cnr <- subset_on_bins(cnr, bins = b2)

    if(!all) {
        cnr[["genes"]] <- cnr$genes[, gg]
        cnr[["gene.index"]] <- cnr$gene.index[cnr$gene.index$hgnc.symbol %in% genes, ]

    }
    
    return(cnr)
} # end subset_on_genes


#' subset based on complete chromosomes
#' @param cnr a cnr bundle
#'
#' @param chrom a list of chromosomes of interst
#' 
#' @keywords internal
#' @noRd
subset_on_chrom <- function(cnr, chrom) {

    b2 <- which(cnr$chromInfo$bin.chrom %in% chrom)
    
    cnr <- subset_on_bins(cnr, bins = b2)
    
    return(cnr)
} # end subset_on_chrom


#' subset based on complete chromosomes
#' @param cnr a cnr bundle
#'
#' @param chrom a chromosome, only one chromosome
#'
#' @param start start coordinate
#'
#' @param end end coordinate
#'
#' @param padding number of bins to add before or after the start and end coordinates
#' 
#' @keywords internal
#' @noRd
subset_on_coordinates <- function(cnr, chrom, start, end, padding = 1) {

    b2 <- which(cnr$chromInfo$bin.chrom %in% chrom &
                cnr$chromInfo$bin.start >= start &
                cnr$chromInfo$bin.end <= end)
    b2 <- c(min(b2) - padding, b2, max(b2) + padding)
    
    cnr <- subset_on_bins(cnr, bins = b2)
    
    return(cnr)
} # end subset_on_coordinates
