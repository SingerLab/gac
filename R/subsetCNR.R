#'  subset a set of bins or genes from a CNR object
#'
#' This function pulls out a section of the genome for further analysis. 
#'
#'  ** still in development **
#'
#' 
#' It currently only works based on bin data, but working to get gene based subset and chromsome position subset.
#' Idealy I'd like to provide coordinate data and bins
#'
#'
#' @param cnr a cnr bundle
#'
#' @param based.on "X", "genes", "position"
#'
#' @param keep string vector with which bins or genes to keep
#'
#' @param ... additional parameters
#' 
#' @return
#'
#' Returns a cnr object with only the desired bins.  This function also subsets the gene index and chromInfo as well
#'
#' @examples
#' data(cnr)
#' data(segCol)
#' 
#' chr19 <- which(as.character(cnr$chromInfo$chrom) == "19")
#' cnr_chr19 <- subsetCNR(cnr, based.on = "X", keep = chr19)
#'
#' sapply(cnr_chr19, dim)
#'
#' HeatmapCNR(cnr_chr19, col = segCol)
#' 
#' @export
subsetCNR <- function(cnr, based.on = "X", keep, ...) {

    message("work in progress, ping me if you really need this")

    kgb <- cnr$gene.index$hgnc.symbol[cnr$gene.index$bin.id %in% keep]
    
    if(based.on == "X") {
        muffin <- cnr$X[keep, ]
        puffin <- cnr$genes[, kgb]
        Ye <- cnr$exprs[, kgb]
        gene.index <- cnr$gene.index[cnr$gene.index$bin.id %in% keep, ]
        chromInfo <- cnr$chromInfo[keep, ]
    }

    cnr[["X"]] <- muffin
    cnr[["Y"]] <- Y
    cnr[["genes"]] <- puffin
    cnr[["exprs"]] <- Ye
    cnr[["gene.index"]] <- gene.index
    cnr[["chromInfo"]] <- chromInfo

    return(cnr)

}
