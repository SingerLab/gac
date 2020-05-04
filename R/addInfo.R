#' Adds info to the chromInfo
#'
#' The function adds information at the bin level.  The gene information is an
#' interpolation of the bin.  The number of independent tests you are performing
#' is the number of common, non-identical segments in the data.   However, to keep
#' the framework simple, it's best to append results from algorithms to the
#' chromInfo. E.g. if .X is 5000 bins, then your genome wide algorithim
#' is 5000 p-values.  Here is where you add it after generating it.
#'
#'
#' @param cnr a cnr bundle
#'
#' @param df a data.frame with the data to incorporate.  Particularly useful for adding p-values, genetic effects, etc
#' 
#' 
#' 
#' @export
addInfo <- function(cnr, df) {
    if(nrow(cnr$chromInfo) == nrow(df)) {
        
        Info <- cbind(cnr$chromInfo, df)
        
    }

    cnr <- list(cnr$X, cnr$Y, cnr$exprs, cnr$qc, Info, cnr$gene.index)
    names(cnr) <- c("X", "Y", "exprs", "qc", "chromInfo", "gene.index")
    
    return(cnr)
}
