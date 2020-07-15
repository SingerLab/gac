#'  subset a set of bins or genes from a CNR object
#'
#' This function pulls out a section of the genome for further analysis. 
#'
#' It currently only works based on bin data, but working to get gene based
#' subset and chromsome position subset.
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
#' chr19 <- rownames(cnr$X)[cnr$chromInfo$chrom == "19"]
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

        ## select hgnc sybmols within bins, and present in `genes`
        kgb <- cnr$gene.index$hgnc.symbol[cnr$gene.index$bin.id %in% keep]
        kgb <- kgb[kgb %in% names(cnr$genes)]

        ## subsetting by bin.id
        muffin <- cnr$X[keep, ]
        chromInfo <- cnr$chromInfo[keep, ]
        gene.index <- cnr$gene.index[cnr$gene.index$bin.id %in% keep, ]
        ## subsetting by hgnc.symbol
        puffin <- cnr$genes[, kgb]
        
        if(!is.null(cnr$exprs)) {
            Ye <- cnr$exprs[, gene.index$hgnc.symbol]
        } else {
            Ye <- NULL
        }
        
        
    } else {
        
        if(based.on == "genes" & all(keep %in% cnr$chromInfo$hgnc.symbol)) {

            kgb <- cnr$chromInfo$bin.id[cnr$chromInfo$hgnc.symbol %in% keep]
            
            muffin <- cnr$X[kgb, ]
            chromInfo <- cnr$chromInfo[kgb, ]
            ## keeping gene index based on bin.id
            gene.index <- cnr$gene.index[cnr$gene.index$bin.id %in% kgb, ]
            puffin <- cnr$genes[, keep]
            
            if(!is.null(cnr$exprs)) {
                Ye <- cnr$exprs[, gene.index$hgnc.symbol]
            } else {
                Ye <- NULL
            }
        }
    }

    cnr[["X"]] <- muffin
    cnr[["genes"]] <- puffin
    cnr[["exprs"]] <- Ye
    cnr[["gene.index"]] <- gene.index
    cnr[["chromInfo"]] <- chromInfo

    return(cnr)

}
