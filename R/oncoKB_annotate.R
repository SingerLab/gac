#' annotate gene index with oncoKB
#' 
#' Cross-references a gene index with OncoKB.  Information about OncoKB is available
#'  at https://www.oncokb.org
#'
#' When using OncoKB, please cite: [Chakravarty et al., JCO PO 2017.](https://ascopubs.org/doi/full/10.1200/PO.17.00011)
#' 
#' @param cnr a cnr bundle
#'
#' @param oncokb an oncokb table
#' 
#' @return
#' Function retunrs an annotated gene.index with OncoKB cancer.gene,oncogene, and
#'  tummor suppressor gene (tsg).  Currently does not bring in the Actionability
#'  OncoKB level annotations.  
#'
#' @examples
#'
#'\dontrun{
#' oncoKB <- read.delim("inst/extdata/oncokb.txt")
#'
#' cnr <- oncoKB_annotate(cnr, oncokb = oncoKB)
#' 
#'}
#' 
#' @keywords internal
#' @noRd
oncoKB_annotate <- function(cnr, oncokb) {

    ngi <- cnr$gene.index
    ngi$oncoKB <- NA
    ngi$oncoKB[ngi$hgnc.symbol %in% oncokb$Hugo.Symbol] <- "cancer.gene"
    ngi$oncoKB[ngi$hgnc.symbol %in% oncokb$Hugo.Symbol[oncokb$OncoKB.Annotated]] <- "OncoKB"
    ngi$oncoKB[ngi$hgnc.symbol %in% oncokb$Hugo.Symbol[oncokb$FOUNDATION.ONE.HEME]] <- "FOUNDATION.ONE.HEME"
    ngi$oncoKB[ngi$hgnc.symbol %in% oncokb$Hugo.Symbol[oncokb$FOUNDATION.ONE]] <- "FOUNDATION.ONE"
    ngi$oncoKB[ngi$hgnc.symbol %in% oncokb$Hugo.Symbol[oncokb$MSK.HEME]] <- "MSK.HEME"
    ngi$oncoKB[ngi$hgnc.symbol %in% oncokb$Hugo.Symbol[oncokb$MSK.IMPACT]] <- "MSK.IMPACT"
    ngi$oncoKB[ngi$hgnc.symbol %in% oncokb$Hugo.Symbol[oncokb$Is.Oncogene]] <- "oncogene"
    ngi$oncoKB[ngi$hgnc.symbol %in% oncokb$Hugo.Symbol[oncokb$Is.Tumor.Suppressor.Gene]] <- "tsg"
    
    cnr[["gene.index"]] <- ngi
    
    return(cnr)
    
}

