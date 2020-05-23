#' Expands the bins to genes
#'
#' This function interpolates the bins into an interpretable gene copy number matrix. All genes within the bin coordinates get assigned the same value. Some genes sit between bins without copy number data and will be marked as NA
#'
#' The gene copy number matrix is core to the cnr object.
#'
#' @param X the X bins data.frame
#'
#' @param gene.index gene index with bin.id and hgnc.symbol
#'
#' @param bin.id name of the bin.id column (the bin.id is the row in the bins data.frame)
#'
#' @param gene.id string noting which gene ID you want to use e.g. hgnc.symbol, ensembl_gene_id
#'
#' 
#'
#' 
#' @export
expand2genes <- function(X, gene.index, bin.id = "bin.id", gene.id = "hgnc.symbol") {
    
    ## expanding segment data to geneCN data
    giu <- gene.index[!duplicated(gene.index[, gene.id]) & gene.index[, gene.id] != "", ]
    geneCN  <- t(X[giu[, bin.id], ])
    colnames(geneCN) <- giu[, gene.id]
    geneCN <- data.frame(geneCN)
    rownames(geneCN) <- colnames(X)
    
    return(geneCN)
    
} ## expand2genes
