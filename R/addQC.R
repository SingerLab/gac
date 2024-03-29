#' addQC
#'
#' Add QC columns to the qc data frame
#'
#' @param cnr a cnr bundle
#'
#' @param df data containg the additional metadata to be incorporated into the qc
#'
#' @param ... additional parameters for sort
#'
#' @return
#'
#' Returns a CNR object with added columns to QC, e.g. mapd, DNA concentratios
#' 
#' @examples
#'
#' data(cnr)
#' 
#' mapd <- data.frame(t(apply(cnr$X, 2, mapd)))
#' mapd <- data.frame(cellID = rownames(mapd), mapd)
#' 
#' cnr <- addQC(cnr, df = mapd, by = "cellID", sort = FALSE)
#' 
#' @export
addQC <- function(cnr, df, ...) {
    
    if(all(df$cellID %in% cnr$qc$cellID)) {
        
        QC <- merge(cnr$qc, df, ...)
        rownames(QC) <- QC$cellID
    }
    
    cnr[["qc"]] <- QC

    cnr <- sync_cnr(cnr)
    
    return(cnr)
}
