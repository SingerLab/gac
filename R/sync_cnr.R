#' sync cnr cells to those in the phenotype annotation
#'
#' @param cnr a cnr bundle
#' @param cell.order a specific order of cells.  `cell.order` must contain all cells.
#'  For subsetting cells use keepCells or excludeCells. Default: NULL which will
#'  syncronize to the in the Y matrix.  
#'
#' @return
#' Function returns a syncronized cnr.
#'
#' @examples
#' data(cnr)
#' names(cnr$X)
#'
#' cnrS <- sync_cnr(cnr)
#' names(cnrS$X)
#' 
#' ordered.cells <- cnr$Y[order(cnr$Y$random1), "cellID"]
#' cnrS <- sync_cnr(cnr, cell.order = ordered.cells)
#' names(cnrS$X)
#' 
#' @importFrom assertthat assert_that
#' @export
sync_cnr <- function(cnr, cell.order = NULL) {
    
    assertthat::assert_that(nrow(cnr$Y) == ncol(cnr$X))
    assertthat::assert_that(nrow(cnr$Y) == nrow(cnr$qc))
    assertthat::assert_that(nrow(cnr$Y) == nrow(cnr$genes))

    if(!is.null(cell.order)) {
        assertthat::assert_that(all(rownames(cnr$Y) %in% cell.order))
        assertthat::assert_that(all(cell.order %in% rownames(cnr$Y)))
        assertthat::assert_that(all(cell.order %in% rownames(cnr$qc)))
        assertthat::assert_that(all(cell.order %in% rownames(cnr$genes)))
        assertthat::assert_that(all(cell.order %in% colnames(cnr$X)))
        
        cnr[["Y"]] <- cnr$Y[cell.order, ]
        rownames(cnr$Y) <- cnr$Y$cellID
        
    } else {

        cell.order <- cnr$Y$cellID
        rownames(cnr$Y) <- cnr$Y$cellID

    }

    if(!is.null(cnr$exprs)) {
        assertthat::assert_that(all(rownames(cnr$exprs) %in% cell.order))
        assertthat::assert_that(all(cell.order %in% rownames(cnr$exprs)))
        cnr[["exprs"]] <- cnr$exprs[cell.order, ]
    }
    
    assertthat::assert_that(all(cell.order %in% colnames(cnr$X)))
    cnr[["X"]] <- cnr$X[, cell.order]
    
    assertthat::assert_that(all(cell.order %in% rownames(cnr$genes)))
    cnr[["genes"]] <- cnr$genes[cell.order, ]

    assertthat::assert_that(all(cell.order %in% cnr$qc$cellID))
    rownames(cnr$qc) <- cnr$qc$cellID
    cnr[["qc"]] <- cnr$qc[cell.order, ]

    cnr[["cells"]] <- cnr$Y$cellID

    return(cnr)
}
