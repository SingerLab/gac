#' pull a set of copy numbers from a cnr into a data.frame
#'
#' @param cnr a cnr bundle
#'
#' @param genes character, genes of interst
#'
#' @return
#' A `data.frame` containing the copy numbers for selected genes.  The data set
#'  is restricted to the genes present in the `genes` matrix.  However, those not
#'  found will be listed out.
#'
#' @examples
#' data(cnr)
#'
#' pull_gene_copy_numbers(cnr, c("MDM2", "CDK4"))
#'
#' genes.of.interest <- list_gene_symbols(cnr, at = "12:58140000:58200000")
#'
#' pull_gene_copy_numbers(cnr, genes.of.interest)
#'
#' @export
pull_gene_copy_numbers <- function(cnr, genes) {
    genesRN <- gsub("-", ".", genes)
    ##
    not.in.set <- genes[! genesRN %in% colnames(cnr[["genes"]])]
    genesRN <- genesRN[genesRN %in% colnames(cnr[["genes"]])]
    ##
    if(length(not.in.set) != 0) { message("Genes ", not.in.set, " were not found") }
    ##
    if(length(genes) == 1) {
        geneCN <- data.frame(cnr[["genes"]][, genesRN])
        names(geneCN) <- genesRN
    } else {
        geneCN <- cnr[["genes"]][, genesRN]
    }
    ##
    return(geneCN)
} # end pull_gene_copy_numbers
