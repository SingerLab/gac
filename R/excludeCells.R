#' exclude cells from CNR bundle
#'
#' Both keepCells and excludeCells perform similar functions; keep cells with positive selection, and excludeCells with negative selection
#' 
#' 
#' @param cnr the cnr bundle
#'
#' @param excl a string vector of cellID to be removed
#'
#'
#' @examples
#'
#' data(cnr)
#'
#' cnr2 <- excludeCells(cnr, excl = cnr$Y$cellID[cnr$qc$ReadsKept < 800000])
#'
#' ## not run: similarly:
#' 
#'
#' cnr3 <- keepCells(cnr, keep = cnr$Y$cellID[cnr$qc$ReadsKept >= 800000])
#'
#' all.equal(cnr2, cnr3)
#'
#' 
#' @export
excludeCells <- function(cnr, excl) {

    if(all(excl %in% colnames(cnr$X))) {
        keep <- colnames(cnr$X)[!colnames(cnr$X) %in% excl]
        X <- cnr$X[,keep ]
        Y <- cnr$Y[keep, ]
        genes <- cnr$genes[keep, ]
        qc <- cnr$qc[keep, ]
        Ye <- cnr$exprs[keep, ]

    } else {

        warning("not all cells in `excl` are in X, removed what was possible")

        keep <- colnames(cnr$X)[!colnames(cnr$X) %in% excl]
        X <- cnr$X[, keep]
        Y <- cnr$Y[keep, ]
        genes <- cnr$genes[keep, ]
        qc <- cnr$qc[keep, ]
        Ye<- cnr$exprs[keep, ]

    }

    cnr <- list(X, genes, Y, qc, Ye, cnr$chromInfo, cnr$gene.index)
    names(cnr) <- c("X", "genes", "Y", "qc", "exprs", "chromInfo", "gene.index")

    return(cnr)
    
}
