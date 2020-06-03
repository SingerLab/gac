#' addQC
#'
#'
#' Add QC columns to the qc data frame
#'
#' @param cnr a cnr bundle
#'
#' @param df data containg the additional metadata to be incorporated into the qc
#'
#' @param by sort by common colum between the cnr$qc and df
#'
#' @param sort logical, weather `merge` sorts the table or not, default is FALSE
#'
#' @param ... additional parameters for sort
#'
#' @return
#'
#' returns a CNR bundle with new QC columns
#' 
#' @examples
#'
#' data(cnr)
#' 
#' mapd <- data.frame(t(apply(cnr$X, 2, mapd)))
#' mapd <- data.frame(cellID = rownames(mapd), mapd)
#' 
#' cnr <- addQC(cnr, df = mapd)
#' 
#' saveRDS(cnr, file = "cnr.rds")
#'
#' 
#' @export
addQC <- function(cnr, df, by = "cellID", sort = FALSE, ...) {
    
    if(all(df$cellID %in% cnr$qc$cellID)) {
        
        QC <- merge(cnr$qc, df, by = by, sort = sort, ...)
        rownames(QC) <- QC$cellID
    }
    
    cnr[["qc"]] <- QC
    
    return(cnr)
}
