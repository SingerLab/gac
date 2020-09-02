#' pull cluster membership for K clusters
#'
#'
#' @param cnr a cnr object
#'
#' @param K number of K clusters to use
#'
#' @export
setKcc <- function(cnr, K = 5) {
    
    cc_membership <- as.data.frame(cnr[["ccp"]][[K]]$consensusClass)
    colnames(cc_membership) <- paste0("ccK", K)

    cnr <- addPheno(cnr, df = cc_membership, by.x = "cellID", by.y = 0,
                    sort = FALSE)
    return(cnr)
}

