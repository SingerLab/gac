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

    binDF <- get_bin_frequencies_cn(cnr)
    geneDF <- get_gene_frequencies_cn(cnr)
    rownames(geneDF) <- gsub("\\.", "-", rownames(geneDF))

    assertthat::assert_that(all(rownames(geneDF) %in%
                                cnr[["gene.index"]][, symbol]))

    cnr <- addInfo(cnr, df = binDF, gdf = geneDF,
                   by.x = symbol, by.y = 0)
    
    ## cnr[["gene.index"]] <- merge(cnr[["gene.index"]], geneDF, by.x = symbol,
    ## by.y = 0, sort = FALSE, all.x = TRUE, ...)
    ## rownames(cnr[["gene.index"]]) <- cnr[["gene.index"]][, symbol]
    
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
get_bin_frequencies_cn <- function(cnr, ...) {

    if(cnr$bulk) {
        X <- cbioportal_states_from_log2_ratio(
            cnr,
            cbioportal =  FALSE, ...)
    }
    if(!cnr$bulk) {
        X <- cnr[["X"]]
    }

    
    AmpCT <- get_amp_counts(X, bulk = cnr$bulk, ...)
    delCT <- get_del_counts(X, bulk = cnr$bulk, ...)
    altCT <- get_alt_counts(X, bulk = cnr$bulk, ...)

    ncells <- ncells(cnr)
    
    altDF <- data.frame(AmpCT, delCT, altCT,
                        AmpFQ = AmpCT/ncells,
                        delFQ = delCT/ncells,
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
get_gene_frequencies_cn <- function(cnr, ...) {

    ## assertthat::assert_that(!cnr$bulk)

    if(cnr$bulk) {
        X <- cbioportal_states_from_log2_ratio(
            cnr, genes =  TRUE,
            cbioportal =  TRUE)
    }
    if(!cnr$bulk) {
        X <- t(cnr[["genes"]])
    }
    
    AmpCT <- get_amp_counts(X, bulk = cnr$bulk, ...)
    delCT <- get_del_counts(X, bulk = cnr$bulk, ...)
    altCT <- get_alt_counts(X, bulk = cnr$bulk, ...)

    nx <- ncells(cnr)
    
    altDF <- data.frame(AmpCT, delCT, altCT,
                        AmpFQ = AmpCT/nx,
                        delFQ = delCT/nx,
                        altFQ = altCT/nx)
    
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
get_amp_counts <- function(X, bulk, base.ploidy = 2,
                           log2.ratio.gain.threshold = 0.4) {
                           
    ampX <- X
    
    if(!bulk) {
        ampCT <- ampX > base.ploidy
    }
    
    if(bulk) {
        ampCT <- ampX >= log2.ratio.gain.threshold
    }
        
    ampCT <- rowSums(ampCT)
    return(ampCT)
    
} ## get_amp_counts


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
get_del_counts <- function(X, bulk,
                           base.ploidy = 2,
                           log2.ratio.loss.threshold = -0.6) {

    delX <- X
    
    if(!bulk) {
        delCT <- delX < 2
    }
    if(bulk) {
        delCT <- delX <= log2.ratio.loss.threshold
    }

    delCT <- rowSums(delCT)
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

    if(bulk) {
        X <- callX(X, ...)
    }
    altX <- binary.X(X, bulk =  bulk, ...)
    
    altCT <- rowSums(altX)
    
    return(altCT)
    
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
callX <- function(X,
                  deletion.threshold =  -1.2,
                  loss.threshold = -0.6,
                  gains.threshold =  0.4,
                  amplification.threshold = 0.8) {
    
    cx <- apply(X, 2, function(i) {
        out <- ifelse(i <= loss.threshold, -1,
               ifelse(i <= deletion.threshold, -2,
               ifelse(i > loss.threshold &
                      i < gains.threshold, 0, 
               ifelse(i >= gains.threshold[1] &
                      i < amplification.threshold[1], 1,
               ifelse(i >= amplification.threshold[1], 2, NA)))))
        })
    
    return(cx)
    
}
