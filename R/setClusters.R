#' setting clusters at optimal height
#'
#' @param cnr a cnr bundle
#'
#' @param tree.height height of the tree (from optClust)
#'
#' @param prefix a prefix to append before the cluster number
#' 
#' @import assertthat
#' @importFrom stats cutree
#' 
#' @export
setClusters <- function(cnr, tree.height, prefix = "C") {

    assertthat::assert_that("cluster" %in% colnames(cnr[["Y"]]) == FALSE,
                            msg = "a cluster column already exists")

    cluster <-  cutree(cnr[["hcdb"]], h = tree.height)

    if(!is.null(prefix)) {
        cluster <- paste0(prefix, cluster)
    }
    
    cnr[["Y"]]$cluster <- cluster
    cnr[["tree.height"]] <- tree.height

    return(cnr)
}
