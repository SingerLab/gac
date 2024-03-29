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

    if(is.null(cnr[["cdb"]])) {
        cnr <- phylo_cnr(cnr)
    }
    
    if(iters <= 201) {
        message("Default value of iters is set to 200. This number iterations only shows general trends.  To identify rare events and off-diagonal events with accuracy, please consider increasing this parameter according to the complexity of your data, and the frequency of events your interested.")
        }
    
    cnr[["ccp"]] <- ConsensusClusterPlus::ConsensusClusterPlus(cnr[["cdb"]],
                                         maxK = maxK, reps = iters,
                                         title = title, innerLinkage = innerLinkage,
                                         finalLinkage = finalLinkage, seed = seed,
                                         verbose = TRUE,
                                         ...)
    return(cnr)
}
