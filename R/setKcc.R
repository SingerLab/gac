#' pull cluster membership for K clusters
#'
#'
#' @param cnr a cnr object
#'
#' @param kCC number of K clusters to use, 
#'  select optimum or maximum k from kSpectral.  Default is NULL
#'  which will pull the optK[kCC] from kStats
#'
#' @param prefix prefix charcter to append to Consensus Clusters
#' 
#' @import assertthat
#' 
#' @return
#'
#' returns cluster membership based on consensus clustering for a
#'  specified kCC
#' 
#' @export
setKcc <- function(cnr, kCC = NULL, prefix = NULL) {

    assertthat::assert_that(any("ConsensusC" %in% colnames(cnr[["Y"]])),
                            msg = "a consensus cluster column already exists")
    
    if(is.null(kCC) & !is.null(cnr[["optK"]]["kCC"])) {
        kCC <- cnr[["optK"]]["kCC"]
        message("Using recommended kCC = ", kCC, " by kSpectral")
    } else {
        assertthat::assert_that(kCC < length(cnr[["ccp"]]))
    }
    
    cc_membership <- as.data.frame(factor(cnr[["ccp"]][[kCC]]$consensusClass))
    colnames(cc_membership) <- "ConsensusC"

  if(!is.null(prefix)) {
        cc_membership[,1] <- paste0(prefix, cc_membership[,1])
    }
    
    cnr <- addPheno(cnr, df = cc_membership, by.x = "cellID", by.y = 0,
                    sort = FALSE)

  cnr[["kCC"]] <- kCC

  return(cnr)
}

