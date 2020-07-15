#' Adds additional phenotypes for a new trait
#'
#' add phenotype columns to the data
#'
#'
#' @param cnr a cnr bundle
#'
#' @param df a data frame with new Y traits
#'
#' @param ... additional arguments by merge
#'
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
addPheno <- function(cnr, df, ...) {

    if(all(df$cellID %in% colnames(cnr$X))) {
        ny <- merge(cnr$Y, df, ...)
        rownames(ny) <- cnr$Y$cellID
    }
    
    cnr[["Y"]] <- ny
    
    return(cnr)
    
}
