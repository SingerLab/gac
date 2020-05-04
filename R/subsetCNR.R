#'  subset a set of bins or genes from a CNR object
#'
#' This function pulls out a section of the genome for further analysis. 
#'
#'  ** still in development **
#'
#' 
#' It currently only works based on bin data, but working to get gene based subset and chromsome position subset.
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
#' chr19 <- which(cnr$chromInfo$chr == "19")
#' cnr_chr19 <- subsetCNR(cnr, based.on = "X", keep = chr19)
#'
#' sapply(cnr_chr19, dim)
#'
#' HeatmapCNR(cnr_chr19, col = segCol)
#' 
#' @export
subsetCNR <- function(cnr, based.on = "X", keep, ...) {

    message("work in progress, ping me if you really need this")

    if(based.on == "X") {
        muffin <- cnr$X[keep, ]
        kgb <- cnr$gene.index$hgnc.symbol[cnr$gene.index$bin.id %in% keep]
        puffin <- cnr$genes[, kgb]
        Y <- cnr$Y
        Ye <- cnr$exprs
        qc <- cnr$qc

        gene.index <- cnr$gene.index[cnr$gene.index$bin.id %in% keep, ]

        chromInfo <- cnr$chromInfo[keep, ]
    }

    cnr <- list(muffin, puffin, Y, Ye, qc, chromInfo, gene.index)
    names(cnr) <- list("X", "genes", "Y", "epxrs", "qc", "chromInfo", "gene.index")

    return(cnr)
    
}
