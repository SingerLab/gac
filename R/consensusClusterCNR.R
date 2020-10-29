#' run consensus clustering
#'
#' Function runs consensus clustering for a pre-determined sequence of K.
#'
#' @param cnr a cnr object
#'
#' @param maxK maximum K, end of K sequence
#'
#' @param reps number of times the clustering runs
#'
#' @param title title for the consensus cluster map
#'
#' @param innerLinkage inner linkage, see ConsensusClusterPlus
#'
#' @param finalLinkage final linkage, see ConsensusClusterPlus
#'
#' @param seed seed of analysis
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
#' cnr <- phyloCNR(cnr)
#'
#' cnr <- consensusClusterCNR(cnr, maxK = 6)
#' 
#' @import ConsensusClusterPlus
#'
#' @export
consensusClusterCNR <- function(cnr, maxK = 40, reps = 150,
                                title = "cnr_ccp",
                                innerLinkage = "ward.D2",
                                finalLinkage = "ward.D2",
                                seed = 2020.0314,
                                ...) {

    assertthat::assert_that(!is.null(cnr[["cdb"]]))

    cnr[["ccp"]] <- ConsensusClusterPlus::ConsensusClusterPlus(cnr[["cdb"]],
                                         maxK = maxK, reps = reps,
                                         title = title, innerLinkage = innerLinkage,
                                         finalLinkage = finalLinkage, seed = seed,
                                         ...)
    return(cnr)
}
