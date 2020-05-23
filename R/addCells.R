#' Add cells to an existing cnr
#'
#' The function deals a little with the n+1 problem.
#'
#' Your data needs to `exact` to the existing cnr for it to be compatible.
#'  This means it needs to have the same bin coordinates and number of X bins, annotation colums, and same qc colums.  It is assumed you generated it with the same pipeline
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
#' @param ... additional parameters if needed
#'
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
#' ## new X w simulated data
#' newX <- data.frame(cbind(rep(c(5,2), c(3000, 2000)),
#'                         rep(c(2,4),  c(3000, 2000))))
#' names(newX) <- paste0("cell", 13:14)
#' head(newX)
#' tail(newX)
#' 
#' ## creating new phenotypes
#' newY <- head(pheno, n = 2)
#' newY$cellID <- paste0("cell", 13:14)
#' rownames(newY) <- newY$cellID
#'
#' newY[, c(6:9)] <- newY[,c(6,9)]+3
#' newY
#' 
#' ## creating new QC
#' newQC <- head(qc, 2)
#' newQC$cellID <- paste0("cell", 13:14)
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
#' saveRDS(cnr, file = "cnr.rds")
#'
#' 
#' @export
addCells <- function(cnr, newX, newY, newqc, newYe = NULL, ...) {

    muffin <- cbind(cnr$X, newX)

    puffin <- rbind(cnr$genes, expand2genes(newX, cnr$gene.index))
    rownames(puffin) <- colnames(muffin)
    
    Y <- rbind(cnr$Y, newY)
    qc <- rbind(cnr$qc, newqc)

    if(is.null(newYe)) {

        Ye <- matrix(NA, nrow = nrow(newY), ncol = nrow(cnr$gene.index),
                     dimnames = list(newY$cellID, cnr$gene.index$hgnc.symbol))
        
    } else {

        Ye <- newYe
        
    }
    
    exprs <- rbind(cnr$exprs, Ye)
    rownames(exprs) <- Y$cellID
    
    cnr <- list(muffin, puffin, Y, exprs, qc, cnr$chromInfo, cnr$gene.index)
    names(cnr) <- c("X", "genes", "Y", "exprs", "qc",  "chromInfo", "gene.index")

    return(cnr)
}
