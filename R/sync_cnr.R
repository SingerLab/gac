#' sync cnr cells to those in the phenotype annotation
#'
#' @param cnr a cnr bundle
#' 
#' @param cell.order a specific order of cells.  `cell.order` must contain all cells.
#' For subsetting cells use keepCells or excludeCells. Default: NULL which will
#'   syncronize to the in the Y matrix.
#'
#' @param full.sync also order chromInfo and gene.index.  default TRUE.
#' Uses bin.chrom, bin.start for ordering bins, and chrom, start for ordering
#'   genes. To use other column names , use full.sync = FALSE, and order bins
#'   and genes separetly.
#' 
#' @param chromosome.order chromosome order to use as levels for arranging chromInfo and gene.index
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
sync_cnr <- function(cnr, cell.order = NULL, full.sync = TRUE,
                     chromosome.order = c(1:22, "X", "Y", "MT")) {
    
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

    if(full.sync) {

        cnr <- order_bins(cnr, chromosome.order = chromosome.order)
        cnr <- order_genes(cnr, chromosome.order = chromosome.order)
        
    }

    return(cnr)
}


#' order genome bins based on chromosome and starting position
#'
#' @param cnr a cnr bundle
#' 
#' @param chromosome.order order for chromosomes, default is 1:22, X, Y, and MT,
#'  corresponding to the human genome
#'
#' @param chrom.column column name for the bin chromosomes. default "bin.chrom"
#'
#' @param start.column column name for bin start. default "bin.start"
#' 
#'
#' @return
#' Function returns an chromInfo ordered by bin chromosomes and start coordinates
#'
#' @examples
#' data(cnr)
#'
#' set.seed(2023)
#' shuffled.bins <- sample(1:nrow(cnr$chromInfo), size = nrow(cnr$chromInfo))
#' cnr$chromInfo <- cnr$chromInfo[shuffled.bins, ]
#' head(cnr$chromInfo)
#' 
#' cnrS <- order_bins(cnr)
#' head(cnrS$chromInfo)
#' 
#' @importFrom assertthat assert_that
#' @export
order_bins  <- function(cnr, chromosome.order = c(1:22, "X", "Y", "MT"),
                       chrom.column = "bin.chrom", start.column = "bin.start") {
    
    cnr$chromInfo[, chrom.column] <- factor(cnr$chromInfo[, chrom.column],
                                         levels = chromosome.order)
    
    cnr$chromInfo[, start.column] <- as.numeric(cnr$chromInfo[, start.column])
    
    nci <- cnr$chromInfo
    nci <- nci[order(nci[, chrom.column], nci[, start.column]), ]
    
    cnr[["chromInfo"]] <- nci
    
    return(cnr)
    
}



#' order gene.index based on chromosome and starting coordinate
#'
#' @param cnr a cnr bundle
#' 
#' @param chromosome.order order for chromosomes, default is 1:22, X, Y, and MT,
#'  corresponding to the human genome
#'
#' @param chrom.column column name for the bin chromosomes. default "chrom"
#'
#' @param start.column column name for bin start. default "start"
#' 
#'
#' @return
#' Function returns an gene.index ordered by chromosomes and start
#'
#' @examples
#' data(cnr)
#'
#' set.seed(2023)
#' shuffled.genes <- sample(1:nrow(cnr$gene.index), size = nrow(cnr$gene.index))
#' cnr$gene.index <- cnr$gene.index[shuffled.genes, ]
#' head(cnr$gene.index)
#' 
#' cnrS <- order_genes(cnr)
#'
#' head(cnrS$gene.index)
#' 
#' @importFrom assertthat assert_that
#' @export
order_genes  <- function(cnr, chromosome.order = c(1:22, "X", "Y", "MT"),
                        chrom.column = "chrom", start.column = "start") {
    
    cnr$gene.index[, chrom.column] <- factor(cnr$gene.index[, chrom.column],
                                         levels = chromosome.order)
    
    cnr$gene.index[, start.column] <- as.numeric(cnr$gene.index[, start.column])
    
    ngi <- cnr$gene.index
    ngi <- ngi[order(ngi[, chrom.column], ngi[, start.column]), ]
    
    cnr[["gene.index"]] <- ngi
    rownames(cnr[["gene.index"]]) <- cnr$gene.index$hgnc.symbol
    
    return(cnr)
    
}

