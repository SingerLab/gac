#' Adds new columns to the chromInfo and gene.index
#'
#' The function adds information at the bin and gene level by running \code{\link{addBinInfo}},
#' and \code{\link{addGeneInfo}}
#'
#' @param cnr a cnr bundle
#'
#' @param df a data.frame with the data to incorporate into chromInfo.  Particularly
#' useful for adding p-values, genetic effects, etc to the bins
#'
#' @param gdf a data.frame with the data to incorporate into gene.index.  Particularly
#' useful for adding p-values, genetic effects, etc to the genes
#' 
#' @param ... other parameters passed to \code{\link{addGeneInfo}}, subsequently to
#' \code{\link[base]{merge}}
#' 
#'
#' @return
#'
#' Returns a CNR object with added columns to the chromInfo and to the gene.index.
#' e.g. p-values for genome wide scans,
#' 
#' 
#' @examples
#'
#' data(cnr)
#'
#' fakeBinPval <- data.frame(pval = runif(5000))
#'
#' fakeGenePval <- data.frame(hgnc.symbol = cnr$gene.index$hgnc.symbol,
#' pval = runif(nrow(cnr$gene.index)))
#'                        
#' cnr <- addInfo(cnr, df = fakeBinPval, gdf = fakeGenePval)
#' 
#' head(cnr$chromInfo)
#' head(cnr$gene.index)
#'
#' ## add info only to bins
#' 
#' fakeBinPval2 <- data.frame(pval = runif(5000))
#' cnr <- addInfo(cnr, df = fakeBinPval2)
#'
#' head(cnr$chromInfo)
#' head(cnr$gene.index)
#' 
#' @export
addInfo <- function(cnr, df, gdf = NULL, ...) {

    out <- addBinInfo(cnr, df = df)

    if(!is.null(gdf))  {
       out <- addGeneInfo(out, df = gdf, ...)
    }
    
    return(out)
}


#' Adds new columns to the chromInfo
#'
#' The function adds information at the bin level.  
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
#' cnr <- addBinInfo(cnr, df = fakePval)
#' 
#' head(cnr$chromInfo)
#'
#' @importFrom assertthat assert_that
#' 
#' @export
addBinInfo <- function(cnr, df) {

    if("bin.id" %in% colnames(df) ) {
        bInfo <- merge(cnr$chromInfo, df, by = "bin.id", sort = FALSE, all.x = TRUE)
        
    } else {

        assertthat::assert_that(nrow(cnr$chromInfo) == nrow(df))

        if(any(names(df) %in% names(cnr$chromInfo))) {
            names(df) <- paste(names(df), "z", sep = ".")
        }
        
        bInfo <- cbind(cnr$chromInfo, df)
    }
    
    cnr[["chromInfo"]] <- bInfo

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
#' @param gene.ids column for gene.ids, must be included in both
#'
#' @param all.index logical, use all of gene.index. When FALSE, return only
#'  the intersect of gene.idex and df.  default TRUE
#'
#' @param ... additional parameters passed to \code{\link[base]{merge}}
#'
#' @return
#' Returns a CNR object with added columns from df to the gene.index.
#'    e.g. p-values for genome wide scans.
#'
#' By default merge occurs using the intersect of column names of
#' the gene.index and df.  To change this, use by.x and by.y.  See
#'  \code{\link[base]{merge}}.  Rownames are added after the merge
#'  from the column specified as gene.ids
#' 
#' @examples
#' data(cnr)
#'
#' fakeGenePval <- data.frame(hgnc.symbol = cnr$gene.index$hgnc.symbol,
#'                        pval = runif(nrow(cnr$gene.index)))
#'
#' cnr <- addGeneInfo(cnr, df = fakeGenePval)
#' 
#' head(cnr$gene.index)
#'
#' @export
addGeneInfo <- function(cnr, df, gene.ids = "hgnc.symbol",
                        sort = FALSE,
                        all.index = TRUE, ...) {
    
    gInfo <- merge(cnr$gene.index, df, sort = sort, all.x = all.index,  ...)
    rownames(gInfo) <- gInfo[, gene.ids]

    cnr[["gene.index"]] <- gInfo
    
    return(cnr)
}
