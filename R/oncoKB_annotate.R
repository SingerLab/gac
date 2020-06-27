#' annotate gene index with oncoKB
#' 
#' Cross-references a gene index with OncoKB.  Information about OncoKB is available at https://www.oncokb.org
#' 
#' @param cnr a cnr bundle
#'
#' @param oncokb an oncokb table
#' 
#' @return
#' Function retunrs an annotated gene.index with OncoKB cancer.gene,oncogene, and tummor suppressor gene (tsg).  Currently does not bring in the Actionability OncoKB level annotations.  
#'
#' @examples
#'
#'\dontrun{
#' oncokb <- read.table("inst/extdata/oncokb.txt")
#'
#' cnr <- oncoKB_annotate(cnr, oncokb)
#' 
#'}
#' 
#' @export
oncoKB_annotate <- function(cnr, oncokb) {

    ngi <- cnr$gene.index
    ngi$oncokb <- "not.oncokb"
    ngi$oncokb[ngi$hgnc.symbol %in% oncokb$Hugo.Symbol] <- "cancer.gene"
    ngi$oncokb[ngi$hgnc.symbol %in% oncokb$Hugo.Symbol[oncokb$Is.Oncogene]] <- "oncogene"
    ngi$oncokb[ngi$hgnc.symbol %in% oncokb$Hugo.Symbol[oncokb$Is.Tumor.Supressor.Gene]] <- "tsg"
    
    cnr[["gene.index"]] <- ngi
    
    return(cnr)
    
}

