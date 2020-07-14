#' Build a CNR bundle (Copy Number, --Rounded)
#'
#' A CNR bundle is a `list` composed of six matrices of class data.frame (mostly). The
#' objective of the CNR is to keep the six matrices syncronized to optimize data
#' management.  It is kept simple in order to easily pull the necesary information
#' to generate a ComplexHetamap.   The functions keepCells, excludeCells, addCells,
#' subsetCNR, manipulate the complete bundle.  Whereas addQC, addPheno, and addInfo
#' manipulate, as the name implies, the QC, Y (phenotype), and chromInfo tables.  
#' 
#'
#' @param X bin or common segment copy number data.  Can be in `numeric`
#' or integer form.  By default it will be rounded by roundCNR.
#'
#' @param Y phenotype and additional cell-level annotation data nrow(Y)
#' needs to equal ncol(X).  There is no check for this yet. requires a
#' `cellID` column
#' 
#' @param qc cell-level quality control metadata e.g. readcount, median
#' bins, mapd, qc.status, or any other data that is technical about the
#' cells. requires a `cellID` column
#'
#' @param exprs slot for expression same-cell (same-sample) gene
#' expression matrix. default is NULL
#' 
#' @param chromInfo bin chromosome and end position in base pairs.
#' Needs to match X
#'
#' @param gene.index a GRanges generated matrix to link bins to genes
#'
#' @param cells a list of cells/samples in the object
#'
#' @examples
#'
#' data(copynumbers)
#' data(pheno)
#' data(qc)
#' data(chromInfo)
#' data(gene.index)
#'
#' cnr <- buildCNR(X = copynumbers, Y = pheno, qc = qc, exprs = NULL,
#' chromInfo = chromInfo, gene.index = gene.index)
#'
#' class(cnr)
#'
#' head(cnr$bins[, 1:5])
#'
#' head(cnr$genes[, 1:5])
#'
#' data(segCol)
#' 
#' HeatmapCNR(cnr)
#' 
#' saveRDS(cnr, file = "cnr.rds")
#'
#' 
#' @export
buildCNR <- function(X, Y, qc, chromInfo, exprs = NULL, gene.index, ...) {
    
    muffin <- roundCNR(X)
    puffin <- data.frame(expand2genes(muffin, gene.index))
    rownames(puffin) <- colnames(muffin)
    rownames(Y) <- Y$cellID
    rownames(qc) <- qc$cellID
    
    if(is.null(exprs)) {

        cnr <-   list(muffin, puffin,  Y,   qc,   chromInfo,   gene.index)
        names(cnr) <- c("X", "genes", "Y", "qc", "chromInfo", "gene.index")

        cnr[["exprs"]] <- NULL
        
    } else {
        Ye <- exprs
        
        cnr <-  list(muffin, puffin,   Y,     Ye,   qc,   chromInfo,   gene.index)
        names(cnr) <- c("X", "genes", "Y", "exprs", "qc", "chromInfo", "gene.index")

    }

    cnr[["cells"]] <- colnames(cnr[["X"]])

    return(cnr)
    
}

