#' Build a CNR bundle (Copy Number, --Rounded)
#'
#' A CNR bundle is a `\code{\link[base]{list}}` composed of six matrices of class data.frame (mostly).
#' The objective of the CNR is to keep single-cell data matrices syncronized 
#' and facilitate interaction with this data.  The functions \code{\link{keepCells}},
#' \code{\link{excludeCells}}, \code{\link{addCells}}, \code{\link{subsetCNR}}, manipulate 
#' the complete bundle. Whereas \code{\link{addQC}}, \code{\link{addPheno}}, and 
#' \code{\link{addInfo}}, as the name implies they add columns to the QC, Y (phenotype), 
#' and chromInfo tables.
#' 
#'
#' @param X bin or common segment copy number data.  Can be in `numeric`
#' or integer form.  By default it will be rounded by \code{\link{roundCNR}}.
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
#' expression matrix. Must have cellID as rownames. default is NULL. 
#' 
#' @param chromInfo bin chromosome and end position in base pairs.
#' Needs to match X
#'
#' @param gene.index a GRanges generated matrix to link bins to genes
#'
#' @param bulk logical specifying if data type is ratio from bulk DNA sequencing
#'   or integer copy number.  If `TRUE` data is an untransformed segment ratio.
#'   If `FALSE` data is integer copy number
#'
#' @param full.sync sync cnr tables to preseve cell order, and sort chromInfo
#'  and gene.index based on chromsome and start position
#'
#' @param chromosome.order chromosome order. Default is a human genome primary
#'  assembly: 1:22, X, Y, and MT.  
#' 
#' @param ... parameters passed to roundCNR
#'
#' @return
#' A CNR bundle composed of genotype, phenotype, and metadata for a single-cell
#' DNA copy number experiment.  See `getting_started.Rmd` vignette for additional
#' details.
#' 
#' @examples
#'
#' ## Example Data
#' ## Copy Number Matrix
#' data(copynumbers)
#' 
#' ## phenotype data
#' data(pheno)
#' 
#' ## Chromosome Information available in Baslan et al. 2012
#' data(chromInfo)
#' 
#' ## gene index see getting started vignette
#' data(grch37.genes.5k)
#' 
#' cnr <- buildCNR(X = copynumbers, Y = pheno, qc = qc, exprs = NULL,
#'                 chromInfo = chromInfo, gene.index = grch37.genes.5k)
#'
#' class(cnr)
#'
#' head(cnr$X[, 1:5])
#'
#' head(cnr$Y[, 1:5])
#'
#' head(cnr$genes[, 1:5])
#'
#' \dontrun{
#'  saveRDS(cnr, file = "cnr.rds")
#' }
#' 
#' @export
buildCNR <- function(X, Y, qc, chromInfo, exprs = NULL, gene.index,
                     bulk = FALSE, full.sync = TRUE,
                     chromosome.order=c(1:22, "X", "Y", "MT"), ...) {

    ## chose if data is ratio or integer CN
    if(bulk) {
        ## for ratio data transfrom to log2
        muffin <- log2(X)
    } else {
        ## round to nearest state (see round CNR)
        muffin <- roundCNR(X, ...)
    }

    assertthat::assert_that(all(c("bin.chrom", "bin.start") %in% names(chromInfo)))
    assertthat::assert_that(all(c("chrom", "start", "hgnc.symbol") %in% names(gene.index)))
    
    ## interpolate to genes
    puffin <- data.frame(expand2genes(muffin, gene.index = gene.index))
    rownames(puffin) <- colnames(muffin)

    ## confirm cellIDs are not duplicated
    assertthat::assert_that(!any(duplicated(Y$cellID)))
    assertthat::assert_that(!any(duplicated(qc$cellID)))

    ## make cellID as key
    rownames(Y) <- Y$cellID
    rownames(qc) <- qc$cellID

    ## build CNR object
    cnr <-   list(muffin, puffin)
    names(cnr) <- c("X", "genes")

    ## coerce order of Y and qc to match. No missing allowed.
    cnr[["Y"]] <- Y[colnames(cnr[["X"]]), ]
    cnr[["qc"]] <- qc[colnames(cnr[["X"]]), ]

    chromInfo <- cbind(bin.id = 1:nrow(chromInfo), chromInfo)
    assertthat::assert_that(nrow(chromInfo) == max(gene.index$bin.id),
                            msg = "number of rows chromInfo does not match gene.index bin.id, please correct the gene.index$bin.id")
    
    cnr[["chromInfo"]] <- chromInfo
    cnr[["gene.index"]] <- gene.index
    rownames(cnr[["gene.index"]]) <- gene.index$hgnc.symbol
    
    ## if expression matrix is available, add it here
    ## must have rownames as cellID/sampleID
    if(is.null(exprs)) {
        cnr[["exprs"]] <- NULL
    } else {
        assertthat::assert_that(all(rownames(exprs) %in% colnames(cnr[["X"]])))
        assertthat::assert_that(all(colnames(exprs) %in%
                                    colnames(cnr[["gene.index"]]$hgnc.symbol)))
        
        cnr[["exprs"]] <- exprs[colnames(cnr[["X"]]), ]
    }
    
    cnr[["cells"]] <- colnames(cnr[["X"]])
    cnr[["bulk"]] <- bulk
    
    cnr <- sync_cnr(cnr, full.sync = full.sync,
                    chromosome.order = chromosome.order)
    
    return(cnr)
    
} ## end buildCNR
