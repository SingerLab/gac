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
#' cnr <- get_alt_frequenciesCNR(cnr)
#'
#' head(cnr$chromInfo)
#' head(cnr$gene.index)
#'
#' @import assertthat
#' 
#' @export
get_alt_frequenciesCNR <- function(cnr, symbol = "hgnc.symbol", ...) {

    binDF <- get_bin_frequencies(cnr)
    geneDF <- get_gene_frequencies(cnr)
    rownames(geneDF) <- gsub("\\.", "-", rownames(geneDF))
    
    cnr <- addInfo(cnr, df = binDF)

    assertthat::assert_that(all(rownames(geneDF) %in% cnr[["gene.index"]][, symbol]))
    cnr[["gene.index"]] <- merge(cnr[["gene.index"]], geneDF, by.x = symbol,
                                 by.y = 0, sort = FALSE, all.x = TRUE, ...)
    
    return(cnr)
}

    

#' estimate bin alteration frequencies
#'
#' @param cnr a cnr bundle
#'
#' @import assertthat
#'
#' @export
get_bin_frequencies <- function(cnr) {

    assertthat::assert_that(!cnr$bulk)
    
    AmpCT <- get_amp_frequencies(t(cnr[["X"]]), bulk = cnr$bulk)
    delCT <- get_del_frequencies(t(cnr[["X"]]), bulk = cnr$bulk)
    altCT <- get_alt_frequencies(t(cnr[["X"]]), bulk = cnr$bulk)

    ncells <- nrow(cnr$X)
    
    altDF <- data.frame(AmpCT, delCT, altCT,
                        AmpFQ = AmpCT/ncells, delFQ = delCT/ncells,
                        altFQ = altCT/ncells)

    return(altDF)
}

#' estimate bin alteration frequencies
#'
#' @param cnr a cnr bundle
#'
#' @export
get_gene_frequencies <- function(cnr) {

    assertthat::assert_that(!cnr$bulk)

    AmpCT <- get_amp_frequencies(cnr[["genes"]], bulk = cnr$bulk)
    delCT <- get_del_frequencies(cnr[["genes"]], bulk = cnr$bulk)
    altCT <- get_alt_frequencies(cnr[["genes"]], bulk = cnr$bulk)

    ncells <- nrow(cnr$X)
    
    altDF <- data.frame(AmpCT, delCT, altCT,
                        AmpFQ = AmpCT/ncells, delFQ = delCT/ncells,
                        altFQ = altCT/ncells)
    
    return(altDF)
}


#' amplification gene frequency
#'
#' @import assertthat
#' 
#' @export
get_amp_frequencies <- function(X, bulk) {

    assertthat::assert_that(!bulk)
    
    ampX <- X
    ampX[X >= 3] <- 1
    ampX[X <= 2] <- 0

    delCT <- colSums(ampX)

    return(delCT)
}


#' deletion frequency
#'
#' @import assertthat
#' 
#' @export
get_del_frequencies <- function(X, bulk) {

    assertthat::assert_that(!bulk)
    
    delX <- X
    delX[X < 2] <- 1
    delX[X >= 2] <- 0
    
    delCT <- colSums(delX)

    return(delCT)
}

#' alteration frequency (both Amp + Del)
#'
#' @import assertthat
#' 
#' @export
get_alt_frequencies <- function(X, bulk) {

    assertthat::assert_that(!bulk)
    
    altX <- binary.X(X)
    
    altCT <- colSums(altX)

    return(altCT)
}
