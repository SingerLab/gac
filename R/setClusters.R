#' setting clusters at optimal height
#'
#' @param cnr a cnr bundle
#'
#' @param tree.height height of the tree (from optClust)
#'
#' @param prefix a prefix to append before the cluster number
#'
#' @return
#'
#' Returns a cnr object with cluster membership based on Bray-Curtis
#'  heirarchical clustering.  Tree height must be specified.
#' 
#' @import assertthat
#' @importFrom stats cutree
#' 
#' @export
setBrayClusters <- function(cnr, tree.height, prefix = "C") {
    
    assertthat::assert_that("BrayC" %in% colnames(cnr[["Y"]]) == FALSE,
                            msg = "a cluster column already exists")
    
    BrayC <-  cutree(cnr[["hcdb"]], h = tree.height)

    if(!is.null(prefix)) {
        BrayC <- paste0(prefix, BrayC)
    }
    
    cnr[["Y"]]$BrayC <- BrayC
    cnr[["tree.height"]] <- tree.height

    return(cnr)
}
