#' annotate gene index with oncoKB
#' 
#' Cross-references a gene index with OncoKB.  Information about OncoKB is available at https://www.oncokb.org
#' 
#'
#' @param cnr a cnr bundle
#'
#' @param oncokb an oncokb table
#'
#' 
#' @return
#' Function retunrs an annotated gene.index with OncoKB cancer.gene,oncogene, and tummor suppressor gene (tsg).  Currently does not bring in the Actionability OncoKB level annotations.  
#'
#'
#' @examples
#'
#'\dontrun{
#' oncokb <- read.table("inst/extdata/oncokb.txt")
#'
#'
#' cnr <- oncoKB_annotate(cnr, oncokb)
#'
#'
#' 
#'}
#' 
#' @export
oncoKB_annotate <- function(cnr, oncokb) {

    gene.index <- cnr$gene.index
    gene.index$oncokb <- "not.oncokb"
    gene.index$oncokb[cnr$gene.index$hgnc.symbol %in% oncokb$Hugo.Symbol] <- "cancer.gene"
    gene.index$oncokb[cnr$gene.index$hgnc.symbol %in% oncokb$Hugo.Symbol[oncokb$Is.Oncogene]] <- "oncogene"
    gene.index$oncokb[cnr$gene.index$hgnc.symbol %in% oncokb$Hugo.Symbol[oncokb$Is.Tumor.Supressor.Gene]] <- "tsg"

    cnr <-     list(X, genes, Y,  exprs, qc, chromInfo,  gene.index)
    names(cnr) <- c("X", "genes", "Y", "exprs", "qc", "chromInfo", "gene.index")

    return(cnr)
    
}


