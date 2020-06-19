#' Build a CNR bundle with bulk DNA data (Copy Number, --Rounded)
#'
#' A CNR bundle is a `list` composed of six matrices of class data.frame (mostly). The
#' objective of the CNR is to keep the six matrices syncronized to optimize data
#' management.  It is kept simple in order to easily pull the necesary information
#' to generate a ComplexHetamap.   The functions keepCells, excludeCells, addCells,
#' subsetCNR, manipulate the complete bundle.  Whereas addQC, addPheno, and addInfo
#' manipulate, as the name implies, the QC, Y (phenotype), and chromInfo tables.  
#'
#'
#' @param X bin or common segment copy number data.  Can be in `numeric` or integer
#' form.  In buildDNA the X matrix is `log2(x+1)` transformed.
#'
#' @param Y phenotype and additional cell-level annotation data nrow(Y) needs to
#' equal ncol(X).  There is no check for this yet. requires a `cellID` column
#' 
#'
#' @param qc cell-level quality control metadata e.g. readcount, median bins, mapd,
#' qc.status, or any other data that is technical about the cells. requires a
#' `cellID` column
#'
#' @param exprs slot for expression same-cell (same-sample) gene expression matrix.
#' default is NULL
#' 
#' @param chromInfo bin chromosome and end position in base pairs. Needs to match X
#'
#'
#' @param gene.index a GRanges generated matrix to link bins to genes 
#'
#' @examples
#' 
#' data(copynumbers)
#' data(pheno)
#' data(qc)
#' data(chromInfo)
#' data(gene.index)
#' 
#' dna <- buildDNA(X = copynumbers, Y = pheno, qc = qc, exprs = NULL,
#' chromInfo = chromInfo, gene.index = gene.index)
#'
#' class(dna)
#'
#' head(dna$X[, 1:5])
#'
#' head(dna$genes[, 1:5])
#'
#' data(lowCol)
#'
#' HeatmapCNR(dna, col = lowCol)
#' 
#' saveRDS(dna, file = "dna.rds")
#'
#' 
#' @export
buildDNA <- function(X, Y, qc, exprs = NULL, chromInfo, gene.index) {
    
    muffin <- log2(X)
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
        
        cnr <-  list(muffin,  puffin,  Y,   Ye,      qc,   chromInfo,   gene.index)
        names(cnr) <- c("X", "genes", "Y", "exprs", "qc", "chromInfo", "gene.index")
    }
    
    return(cnr)
    
}

