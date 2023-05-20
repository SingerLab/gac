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
#' @param df a data.frame with the data to incorporate.  Particularly
#' useful for adding p-values, genetic effects, etc to the bins
#'
#' @return
#'
#' Returns a CNR object with added columns to the chromInfo. e.g. p-values for genome wide scans
#' 
#' @examples
#'
#' data(cnr)
#'
#' fakePval <- data.frame(pval = runif(5000))
#'
#' cnr <- addInfo(cnr, df = fakePval)
#' 
#' head(cnr$chromInfo)
#'
#' @importFrom assertthat assert_that
#' 
#' @export
addInfo <- function(cnr, df) {
    
    if("bin.id" %in% colnames(df) ) {
        Info <- merge(cnr$chromInfo, df, by = "bin.id", sort = FALSE)
    } else {
        assertthat::assert_that(nrow(cnr$chromInfo) == nrow(df))
        Info <- cbind(cnr$chromInfo, df)
    }
    
    cnr[["chromInfo"]] <- Info
    
    return(cnr)
}


#' Adds Gene Information to the gene.index
#'
#' The function adds information to the gene.index.#'
#'
#' @param cnr a cnr bundle
#'
#' @param df a data.frame with the data to incorporate.  Particularly
#' useful for gene annotations, filters e.g. OncoKB, or results from
#' analyses e.g. p-values, genetic effects, etc to the genes
#'
#' @param sort wether to sort the ouput object, default is FALSE
#'
#' @param ... additional parameters passed to merge
#'
#' @return
#' Returns a CNR object with added columns from df to the gene.index.
#'   e.g. p-values for genome wide scans
#' 
#' @examples
#'#' data(cnr)
#'
#' fakePval <- data.frame(pval = runif(5000))
#'
#' cnr <- addGeneInfo(cnr, df = fakePval)
#' 
#' head(cnr$gene.index)
#'
#' @export
addGeneInfo <- function(cnr, df, sort = FALSE, ...) {
    
    gInfo <- merge(cnr$gene.index, df, sort = sort,  ...)
    
    cnr[["gene.index"]] <- gInfo
    
    return(cnr)
}
