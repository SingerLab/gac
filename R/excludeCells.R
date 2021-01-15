#' exclude cells from CNR bundle
#'
#' Both keepCells and excludeCells perform similar functions; keep cells with positive selection, and excludeCells with negative selection
#' 
#' @param cnr the cnr bundle
#'
#' @param excl a string vector of cellID to be removed
#'
#' @return
#'
#' Returns a cnr object after removal of cells
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
        Ye <- cnr$exprs[keep, ]
        
    }
    
    cnr[["X"]] <- X
    cnr[["genes"]] <- genes
    cnr[["Y"]] <- Y
    cnr[["qc"]] <- qc
    cnr[["exprs"]] <- Ye
    cnr[["cells"]] <- keep

    if(!is.null(cnr[["vdj.cells"]])) {
        vk <- cnr[["vdj.cells"]] %in% keep
        cnr[["vdj.cells"]] <- cnr[["vdj.cells"]][vk]
    }

    return(cnr)
    
}
