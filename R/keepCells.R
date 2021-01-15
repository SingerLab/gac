#' exclude cells from CNR bundle
#'
#' Both keepCells and excludeCells perform similar functions; keep cells with positive selection, and excludeCells with negative selection
#' 
#' 
#' @param cnr the cnr bundle
#'
#' @param keep a string vector of cellID to keep
#'
#' @return
#'
#' Returns a CNR bundle with only the cells we want to keep
#'
#' @examples
#'
#' data(cnr)
#'
#' cnr2 <- excludeCells(cnr, excl = cnr$Y$cellID[cnr$qc$ReadsKept < 800000])
#'
#' cnr3 <- keepCells(cnr, keep = cnr$Y$cellID[cnr$qc$ReadsKept >= 800000])
#'
#' all.equal(cnr2, cnr3)
#'
#' 
#' @export
keepCells <- function(cnr, keep) {

    if(all(keep %in% colnames(cnr$X))) {
        X <- cnr$X[, keep]
        Y <- cnr$Y[keep,]
        genes <- cnr$genes[keep,]
        qc <- cnr$qc[keep, ]
        Ye <- cnr$exprs[keep,]

    } else {
        warning("Not all cells are present in the `X` matrix")
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
