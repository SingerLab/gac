#' Compute pin matrix and Fisher distances
#'
#' Builds a breakpoint-based pin matrix from integer copy number profiles,
#' runs a Monte Carlo Fisher test to compute pairwise cell similarity, and
#' stores FDR estimates and an annotated hierarchical tree.  Results are
#' kept in \code{cnr$pins} and the Fisher p-value distance (range 0--1) in
#' \code{cnr$dists$fisher}, ready for combination with Bray-Curtis via
#' \code{\link{phylo_cnr}} (\code{dist.method = "combined"}).
#'
#' @param cnr a cnr bundle (integer copy number, \code{bulk = FALSE})
#'
#' @param cytobands UCSC cytoBand data.frame (columns: chrom, chromStart,
#'   chromEnd, name, gieStain).  Use the bundled \code{hg19_cytoBand} for
#'   human hg19 data.
#'
#' @param centromere cytoband name patterns that define the centromere region
#'   to exclude.  Default \code{c("p11", "q11")}.
#'
#' @param nsim number of Monte Carlo simulations for the null distribution.
#'   Default 200; increase for rare events.
#'
#' @param hc.method linkage method for the annotated hierarchical tree stored
#'   in \code{cnr$pins$hc}.  Default \code{"average"}.
#'
#' @param ... additional arguments passed to \code{.sim_fisher_wrapper}
#'   (e.g. \code{njobs}, \code{nsweep}, \code{seedme}).
#'
#' @return
#' The cnr bundle with new slots:
#' \describe{
#'   \item{\code{cnr$pins}}{Named list: \code{pinmat} (pins x cells binary
#'     matrix), \code{pins} (pin metadata), \code{cells}, \code{ploidies},
#'     \code{centroareas}, \code{centrobins}, \code{fisher_raw} (raw
#'     simulation output), \code{fdr} (log10 FDR matrix), \code{hc}
#'     (annotated hclust tree with FDR and sharing per node).}
#'   \item{\code{cnr$dists$fisher}}{Fisher p-value distance matrix of class
#'     \code{dist}, values in [0, 1]: 0 = maximally similar, 1 = no shared
#'     breakpoint signal.}
#' }
#'
#' @examples
#' \dontrun{
#' data(cnr, hg19_cytoBand)
#' cnr <- sim_fisherCNR(cnr, cytobands = hg19_cytoBand)
#' cnr <- phylo_cnr(cnr, dist.method = "combined")
#' }
#'
#' @importFrom assertthat assert_that
#'
#' @export
sim_fisherCNR <- function(cnr, cytobands, centromere = c("p11", "q11"),
                          nsim = 200, hc.method = "average", ...) {

    assertthat::assert_that(!cnr$bulk,
        msg = "sim_fisherCNR requires integer copy number (bulk = FALSE)")

    centroareas <- .pin_centroareas(cytobands, centromere = centromere)
    centrobins  <- .pin_regions2bins(cnr$chromInfo, centroareas)

    pinX <- cbind(cnr$chromInfo[, c("chrom", "chrompos", "abspos")], cnr$X)
    pins <- .calc_pinmat(cnr$chromInfo, pinX, dropareas = centroareas)
    pins$centroareas <- centroareas
    pins$centrobins  <- centrobins

    pins$fisher_raw <- .sim_fisher_wrapper(pins$pinmat, pins$pins, nsim = nsim, ...)
    if (is.null(pins$fisher_raw)) {
        warning("sim_fisherCNR: fewer than 2 pin sign groups; Fisher test skipped.")
        cnr$pins <- pins
        return(cnr)
    }

    pins$fdr <- .fisher_fdr(pins$fisher_raw$true, pins$fisher_raw$sim, cnr$cells)

    ## annotated tree uses log10(p) internally (more negative = more similar)
    fisher_log_dist <- .cell2cell_matrix(log10(pins$fisher_raw$true), cnr$cells)
    pins$hc <- .hclust_tree(pins$pinmat, pins$fdr, fisher_log_dist,
                            hcmethod = hc.method)

    cnr$pins <- pins

    ## p-values directly as [0,1] distance for combination with Bray-Curtis
    cnr$dists$fisher <- stats::as.dist(
        .cell2cell_matrix(pins$fisher_raw$true, cnr$cells))

    return(cnr)
}
