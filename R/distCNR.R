#' calculating cell-to-cell distances for clustering
#'
#' @param cnr a cnr object
#'
#' @param method method for calculating distances, defaults to Bray-Curtis
#' dissimilarity
#'
#' @param ... other parameters passed to vegan::vegdist
#'
#' @import vegan
#' 
#' @return
#' Returns a cell-to-cell distance matrix of class `dist`
#' 
#' @export
distCNR <- function(cnr, method = "bray", ...) {
    
    if(cnr$bulk) {
        
        ## Bray-Curtis dissimilarity doesn't work well on log2 ratio,
        ## transforming back to ratio
        cdb <- vegan::vegdist(t(2^cnr[["X"]]), method = method, ...)

    } else {

        cdb <- vegan::vegdist(t(cnr[["X"]]), method = method, ...)

        return(cdb)

    }
} ## distCNR


#' heirarchical clustering
#'
#' @param cnr a cnr object
#'
#' @param method method for heirarchical clustering, defaults to "ward.D2"
#'
#' @param ... other parameters passed to hclust
#'
#' @return
#' Returns a heirarchical clustering object
#'
#' @importFrom stats hclust
#' 
#' @export
hclustCNR <- function(cnr, method = "ward.D2", ...) {

    hcdb <- hclust(cnr[["cdb"]], method = method, ...)
    
    return(hcdb)

} ## hclustCNR

