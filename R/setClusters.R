#' setting clusters at optimal height
#'
#' @param cnr a cnr bundle
#'
#' @param tree.height height of the tree (from optClust)
#'
#' @import assertthat
#' 
#' @export
setClusters <- function(cnr, tree.height, prefix = "C") {
    
    cluster <-  cutree(cnr[["hcdb"]], h = tree.height)
    cluster <- paste0(prefix, cluster)

    assertthat::assert_that("cluster" %in% colnames(cnr[["Y"]]) == FALSE)
    
    cnr[["Y"]]$cluster <- cluster

    return(cnr)
}
