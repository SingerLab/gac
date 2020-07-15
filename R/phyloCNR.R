#' calculating distances, heirarchical clustering, and `phylo` object
#'
#' @param cnr a cnr bundle
#'
#' @param dist.method method for calculating cell-to-cell distance
#' (see vegan::vegdist)
#'
#' @param hclust.method method for heirarchical clustering (see hclust)
#'
#' @param ... other parameters
#'
#' @import ape
#'
#' @export
phyloCNR <- function(cnr, dist.method = "bray", hclust.method = "ward.D2", ...) {

    cnr[["cdb"]] <- distCNR(cnr, method = dist.method, ...)

    cnr[["hcdb"]] <- hclustCNR(cnr, method = hclust.method, ...)

    cnr[["phylo"]] <- ape::as.phylo(cnr[["hcdb"]])

    return(cnr)

}

