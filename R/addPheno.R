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
#' rand3 <- cbind(cnr$Y, rnorm(nrow(cnr$Y), mean = 2, sd = 1))
#' 
#' cnr <- addPheno(cnr, df = rand3)
#' 
#' @export
addPheno <- function(cnr, df, ...) {

    if(all(df$cellID %in% colnames(cnr$X))) {

        ny <- merge(cnr$Y, df, by = "cellID", ...)
        rownames(ny) <- cnr$Y$cellID
    }
    
    cnr <- list(cnr$X, ny, cnr$exprs, cnr$qc, cnr$chromInfo, cnr$gene.index)
    names(cnr) <- c("X", "Y", "exprs", "qc", "chromInfo", "gene.index")
    return(cnr)
    
}
