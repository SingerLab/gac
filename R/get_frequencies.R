#' estimate bin and gene alteration frequencies
#'
#' Wrapper function to calculate alteration frequencies (amplifications, deletions, and both combined) at
#' the bin level, and gene level.
#' 
#'
#' @param cnr a cnr bundle
#'
#' @param symbol gene symbol used on cnr[["genes"]] as column names, e.g. HGNC symbol or ensembl ID.
#' defaults to hgnc.symbol on package gene.index
#'
#' @param ... other parameters passed to merge
#'
#' @return
#'
#' Calls on get_bin_frequencies and get_gene_frequencies, respectively, to estimate amplification frequencies,
#' deletion frequencies, and the total of both alterations (amp + del combined) at every bin, and gene.
#'
#' It then merges the bin frequencies to the chromInfo, and gene frequencies to the gene.index.
#'
#' Currently only tested on integer copy number, where amplifications are hard-coded as three (3) or more copies,
#' and deletions as one (1) or zero (0).
#'
#' @examples
#'
#' data(cnr)
#'
#' head(cnr$gene.index)
#' head(cnr$chromInfo)
#'
#' cnr <- get_alteration_frequencies(cnr)
#'
#' head(cnr$chromInfo)
#' head(cnr$gene.index)
#'
#' @importFrom assertthat assert_that
#' 
#' @export
get_alteration_frequencies <- function(cnr, symbol = "hgnc.symbol", ...) {

    binDF <- get_bin_frequencies(cnr)
    geneDF <- get_gene_frequencies(cnr)
    rownames(geneDF) <- gsub("\\.", "-", rownames(geneDF))

    assertthat::assert_that(all(rownames(geneDF) %in%
                                cnr[["gene.index"]][, symbol]))

    cnr <- addInfo(cnr, df = binDF)
    cnr[["gene.index"]] <- merge(cnr[["gene.index"]], geneDF, by.x = symbol,
                                 by.y = 0, sort = FALSE, all.x = TRUE, ...)
    rownames(cnr[["gene.index"]]) <- cnr[["gene.index"]][, symbol]
    
    return(cnr)
}

    

#' estimate bin alteration frequencies
#'
#' @param cnr a cnr bundle
#'
#' @param ... additional parameters passed to get_amp_counts, get_del_counts, or
#' get_alt_counts
#'
#' @importFrom assertthat assert_that
#' 
#' @keywords internal
#' @noRd
get_bin_frequencies <- function(cnr, ...) {

    assertthat::assert_that(!cnr$bulk)
    
    AmpCT <- get_amp_counts(t(cnr[["X"]]), bulk = cnr$bulk, ...)
    delCT <- get_del_counts(t(cnr[["X"]]), bulk = cnr$bulk, ...)
    altCT <- get_alt_counts(t(cnr[["X"]]), bulk = cnr$bulk, ...)

    ncells <- ncol(cnr$X)
    
    altDF <- data.frame(AmpCT, delCT, altCT,
                        AmpFQ = AmpCT/ncells, delFQ = delCT/ncells,
                        altFQ = altCT/ncells)

    return(altDF)
}

#' estimate bin alteration frequencies
#'
#' @param cnr a cnr bundle
#'
#' @param ... additional parameters passed to get_amp_counts, get_del_counts, or
#' get_alt_counts
#' 
#' @importFrom assertthat assert_that
#' 
#' @keywords internal
#' @noRd
get_gene_frequencies <- function(cnr, ...) {

    assertthat::assert_that(!cnr$bulk)

    AmpCT <- get_amp_counts(cnr[["genes"]], bulk = cnr$bulk, ...)
    delCT <- get_del_counts(cnr[["genes"]], bulk = cnr$bulk, ...)
    altCT <- get_alt_counts(cnr[["genes"]], bulk = cnr$bulk, ...)

    ncells <- ncol(cnr$X)
    
    altDF <- data.frame(AmpCT, delCT, altCT,
                        AmpFQ = AmpCT/ncells, delFQ = delCT/ncells,
                        altFQ = altCT/ncells)
    
    return(altDF)
}


#' amplification gene counts
#'
#' @param X a matrix of integer copy number `cnr$X`
#'
#' @param bulk logical, weather the cnr is from bulk DNA, default is FALSE
#'
#' @param base.ploidy tumor expected ploidy
#' 
#' @importFrom assertthat assert_that
#'
#' @keywords internal
#' @noRd
get_amp_counts <- function(X, bulk, base.ploidy = 2) {

    assertthat::assert_that(!bulk)
    
    ampX <- X
    ampX[X >= (base.ploidy+1)] <- 1
    ampX[X <= base.ploidy] <- 0

    delCT <- colSums(ampX)

    return(delCT)
}


#' deletion counts
#'
#' @param X a matrix of integer copy number `cnr$X`
#'
#' @param bulk logical, weather the cnr is from bulk DNA, default is FALSE
#'
#' @param base.ploidy tumor expected ploidy
#'
#' @importFrom assertthat assert_that
#'
#' @keywords internal
#' @noRd
get_del_counts <- function(X, bulk, base.ploidy = 2) {

    assertthat::assert_that(!bulk)
    
    delX <- X
    delX[X < base.ploidy] <- 1
    delX[X >= base.ploidy] <- 0
    
    delCT <- colSums(delX)

    return(delCT)
}

#' alteration counts (both Amp + Del)
#'
#' @param X a matrix of integer copy number `cnr$X`
#'
#' @param bulk logical, weather the cnr is from bulk DNA, default is FALSE
#'
#' @param ... additional parameters passed to get_amp_counts, get_del_counts, or
#' get_alt_counts
#'
#' @importFrom assertthat assert_that
#' 
#' @keywords internal
#' @noRd
get_alt_counts <- function(X, bulk, ...) {

    assertthat::assert_that(!bulk)
    
    altX <- binary.X(X, ...)
    
    altCT <- colSums(altX)

    return(altCT)
}
