#' binaryDDRC
#'
#' @param ddrc a DDRC.df or DDRC.g matrix belonging to a cnr object
#'
#' @param base.ploidy expected cell ploidy, e.g. 2N = 2, 4N = 4.
#'
#' @return
#' 
#' returns a binary/incidence matrix for a DDRC object.  A value of 1 is
#'  returned for all values not equal to the cell ploidy.  Currently being
#' used to create an infScite table using clones. 
#'
#' @examples
#' 
#' data(cnr)
#' 
#' noisy.cells <- cnr$qc$cellID[cnr$qc$qc.status == "FAIL"]
#'
#' ## reduced pipeline to genrate DDRC clone profiles
#' cnr <- excludeCells(cnr, excl = noisy.cells)
#' cnr <- phyloCNR(root.cell = "cell0")
#' cnr <- setBrayClusters()
#' cnr <- consensusClusterCNR(iters = 20, maxK = 40)
#' cnr <- doKSpectral()
#' cnr <- setKcc()
#' cnr <- cluster_heterogeneity(by = "category1",
#'           cluster_column = "ConsensusC")
#' cnr <- get_cluster_profiles()
#' 
#' binary.ddrc <- binaryDDRC(cnr$DDRC.df)
#' head(binary.ddrc)
#' 
#' @export
binaryDDRC <- function(ddrc, base.ploidy = 2) {
    z = matrix(0, nrow = nrow(ddrc), ncol = ncol(ddrc),
               dimnames = list(rownames(ddrc), colnames(ddrc)))
    z[ddrc != base.ploidy] <- 1
    z <- as.data.frame(z)
    z
} ## binaryDDRC

