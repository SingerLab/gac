#' run consensus clustering
#'
#' Function runs consensus clustering for a pre-determined sequence of K.
#'
#' @param cnr a cnr object
#'
#' @param maxK maximum K, end of K sequence
#'
#' @param iters number of times the clustering runs
#'
#' @param title title for the consensus cluster map
#'
#' @param innerLinkage inner linkage, see ConsensusClusterPlus
#'
#' @param finalLinkage final linkage, see ConsensusClusterPlus
#'
#' @param seed seed of analysis
#'
#' @param verbose print out progress
#'
#' @param ... additional parameters for ConsensusClusterPlus e.g. plot = "png"
#'
#' 
#' @return
#'
#' A cnr object that contains the output of ConsensusClusterPlus for the
#'  specified Ks.
#'
#' \itemize{
#'   \item ccp ConsensusClusterPlus object
#' }
#' 
#' 
#' @examples
#' data(cnr)
#'
#' cnr <- phylo_cnr(cnr)
#'
#' cnr <- run_consensus_clustering(cnr, maxK = 6)
#' 
#' @importFrom ConsensusClusterPlus ConsensusClusterPlus
#' @importFrom assertthat assert_that
#'
#' @export
run_consensus_clustering <- function(cnr, maxK = 40, iters = 200,
                                title = "cnr_ccp",
                                innerLinkage = "ward.D2",
                                finalLinkage = "ward.D2",
                                seed = 2020.0314,
                                verbose = TRUE,
                                ...) {

    if (is.null(cnr$dists$bray))
        cnr <- phylo_cnr(cnr)

    ## Use combined distance when available, otherwise Bray-Curtis
    d <- if (!is.null(cnr$dists$combined)) cnr$dists$combined else cnr$dists$bray

    if (iters <= 201)
        message("Default value of iters is set to 200. This number of iterations ",
                "only shows general trends. Consider increasing for rare events.")

    cnr[["ccp"]] <- ConsensusClusterPlus::ConsensusClusterPlus(
        d, maxK = maxK, reps = iters,
        title = title, innerLinkage = innerLinkage,
        finalLinkage = finalLinkage, seed = seed,
        verbose = verbose, ...)

    return(cnr)
}
