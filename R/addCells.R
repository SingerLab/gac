#' Add cells to an existing cnr
#'
#' The function deals a little with the n+1 problem.
#'
#' Your data needs to be `exact` to the existing cnr for it to be compatible.
#'  This means it needs to have the same bin coordinates and number of X
#'  bins, annotation colums, and same qc colums.  It is assumed you generated
#'  it with the same pipeline.
#'
#' @param cnr a cnr bundle
#'
#' @param newX bin data for the new cells to be incorporated
#'
#' @param newY phenotype data for the new cells
#'
#' @param newqc qc metadata for the new cells
#'
#' @param newYe exprssion data for new cells
#'
#' @param do.clean weather to remove additional data stored in the cnr, default is TRUE
#' 
#' @param ... additional parameters if needed
#'
#' @return
#'
#' Returns a CNR object after adding new cells.  If `do.clean = TRUE` only the original
#' eight data items will be returned.
#'
#' \itemize{
#'   \item X copy number data
#'   \item Y phenotype data
#'   \item exprs expression data
#'   \item qc quality control data
#'   \item chromInfo bin to chromosome information
#'   \item gene.index gene to bin and chromosome information
#'   \item cells list of cells
#'   \item bulk bulk DNA or single-cell copy number. If TRUE, data is bulk DNA log2 ratio
#'   \item ... additional analyzed results
#' }
#'
#' @examples
#'
#' library(gac)
#' data(cnr)
#' sapply(cnr, dim)
#'
#' data(pheno)
#' data(qc)
#'
#' ## create unique cell names
#' n.cells <- nrow(cnr$Y)
#' new.cells <- paste0("cell", n.cells + 1:2)
#' 
#' ## new X w simulated data
#' newX <- data.frame(cbind(rep(c(5,2), c(3000, 2000)),
#'                         rep(c(2,4),  c(3000, 2000))))
#' names(newX) <- new.cells
#' head(newX)
#' tail(newX)
#'
#' ## creating new phenotypes
#' newY <- head(pheno, n = 2)
#' newY$cellID <- new.cells
#' rownames(newY) <- newY$cellID
#' newY[, c(6:9)] <- newY[,c(6,9)]+3
#'
#' newY
#' 
#' ## creating new QC
#' newQC <- head(qc, 2)
#' newQC$cellID <- new.cells
#' rownames(newQC) <- newQC$cellID
#'
#' newQC[,2:4] <- newQC[,2:4] + 2
#' newQC
#' 
#' ## add cells
#' cnr <- addCells(cnr, newX = newX, newY = newY, newqc = newQC)
#' sapply(cnr, dim)
#' 
#'
#' @export
addCells <- function(cnr, newX, newY, newqc, newYe = NULL, do.clean = TRUE, ...) {

    muffin <- cbind(cnr$X, newX)

    puffin <- expand2genes(muffin, cnr$gene.index)
    rownames(puffin) <- colnames(muffin)
    
    Y <- rbind(cnr$Y, newY)
    qc <- rbind(cnr$qc, newqc)

    newCells <- names(muffin)
    
    if(!is.null(newYe)) {
        Ye <- rbind(cnr$exprs, newYe)
        rownames(Ye) <- colnames(muffin)
    } else {
        Ye <- NULL
    }
        
    
    if(do.clean) {
        ## re-create cnr
        cnr <- list(muffin, puffin, Y, qc, Ye,
                    cnr$chromInfo, cnr$gene.index, newCells, cnr$bulk)
        names(cnr) <- c("X", "genes", "Y", "qc", "exprs",
                        "chromInfo", "gene.index", "cells", "bulk")
        
    } else {

        cnr[["X"]] <- muffin
        cnr[["genes"]] <- puffin
        cnr[["Y"]] <- Y
        cnr[["exprs"]] <- Ye
        cnr[["qc"]] <- qc
        cnr[["cells"]] <- newCells
    }

    cnr <- sync_cnr(cnr)

    return(cnr)
}
