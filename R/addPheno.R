#' Adds additional phenotypes for a new trait
#'
#' add phenotype columns to the data
#'
#'
#' @param cnr a cnr bundle
#'
#' @param df a data frame with new Y traits
#'
#' @param sort logical, sort data.frame after merge \code{\link[base]{merge}}
#'
#' @param ... additional arguments by merge
#'
#' @return
#'
#' Returns a CNR object with added columns to Y, e.g. additional cell
#'  phenotypes, histologies, etc.
#' 
#' @examples
#'
#' data(cnr)
#' 
#' rand3 <- data.frame(cellID = cnr$Y$cellID,
#'                     rand3 = rnorm(nrow(cnr$Y), mean = 2, sd = 1))
#'                     
#' cnr <- addPheno(cnr, df = rand3, by = "cellID")
#' 
#' @export
addPheno <- function(cnr, df,
                     sort = FALSE, ...) {

    if(all(df$cellID %in% colnames(cnr$X))) {
        ny <- merge(cnr$Y, df, sort = sort, ...)
        rownames(ny) <- ny$cellID
    }
    
    cnr[["Y"]] <- ny

    cnr <- sync_cnr(cnr)
    
    return(cnr)
    
} ## end addPheno
