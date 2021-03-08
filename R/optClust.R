#' optimizing clustering for single-cell copy number
#'
#' The optClust is a simple attempt to set a threshold on where to cut the an
#' hclust tree.  The aim is to maximize the number of cluster with multiple-cells
#' while minimizing the number of one-cell clusters.
#'
#' @return
#'
#' Returns a matrix where rows are the heights in the opt.range, and three columns.
#'
#' One-cell specifies the number of clusters with only one cell, 
#' Multi-cell specifyies the number of clusters with multiple cells (>=2), and
#' %CMC which is the percentage of cells in multi-cell clusters
#'
#' A good threshold is one that minimizes the one-cell cluster, and maximizes
#' the %CMC
#' 
#'
#' @param cnr a cnr bundle
#'
#' @param opt.range  range of tree heights that need to be optimized
#'
#' @examples
#'
#' data(cnr)
#'
#' cnr <- phyloCNR(cnr)
#'
#' optClust(cnr, opt.range = seq(0, 0.3, by = 0.05))
#'
#' #for micro-optimization
#'
#' optClust(cnr, opt.range = seq(0, 0.2, by = 0.001))
#'
#' @export
optClust <- function(cnr, opt.range = seq(0.005, 0.6, by = 0.005)) {

    clL <- sapply(opt.range, function(h) {
        ctbl <- sort(table(cutree(cnr[["hcdb"]], h = h)),
                     decreasing = TRUE)
    })
    names(clL) <- opt.range

    omt <- matrix(t(sapply(clL, function(i) {
        cbind(sum(i == 1), sum(i != 1), sum(i != 1)/length(i))
    })),
    ncol = 3,
    dimnames = list(opt.range, c("One-cell", "Multi-cell", "%CMC")))
    
    return(omt)
    
} # optClust


#' Minimum tree.height at the intersect of one-cell and multi-cell clusters
#'
#' @param mocp output from \link{optClust}
#'
#' @export
minimum.intersect <- function(mocp) {

    mi <- as.numeric(names(which(mocp[,1] < mocp[,2])[1]))
    return(mi)
    
} #end minimum.intersect
