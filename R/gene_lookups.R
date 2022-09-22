#' finds genes within one chromosome interval
#' @param cnr a cnr bundle
#'
#' @param at character, coordinate string
#'
#' @param show.columns subset of columns from gene.index to output.
#'  default is NULL, which shows all columns
#'
#' @param identifier name of column with gene identifiers
#'
#' @return
#' Returns the subset of the gene.index within the coordinates specified.
#'  By default all columns are returned, with the option to select a
#'  subset using \code{show.columns}.
#'
#' @examples
#'
#' data(cnr)
#'
#' find_genes(cnr, at = "12:58,000,000-59,000,000")
#'
#' find_genes(cnr, at = "12:58,000,000-59,000,000",
#'            show.columns = c("hgnc.symbol", "gene.type"))
#' 
#' 
#' @export
find_genes <- function(cnr, at,
                       show.columns = NULL,
                       identifier = "hgnc.symbol") {

    gene.list <- list_genes(cnr = cnr, at = at, identifier = identifier)
    idx <- cnr$gene.index[, identifier] %in% gene.list
    
    if(is.null(show.columns)) {
        out <- cnr$gene.index[idx, ]
    } else {
        assertthat::assert_that(all(show.columns %in% names(cnr$gene.index)))
        out <- cnr$gene.index[idx, show.columns]
    }
    
    return(out)
}

#' convert ucsc style coordinates to ensembl
#' @param x character, coordinate string such as "chr0:000,000,000-000,000,001"
#'   or 0:000,000,000-000,000,000
#'
#' @return
#' Returns a the same character string with the `chr`, and commas `,` removed, 
#'  and all dividers as `:`
#'
#' @examples \dontrun{
#'
#' convert_coord("chr12:58,000,000-59,000,000")
#'
#' }
#'
#' @keywords internal
#' @noRd
convert_coord <- function(x) {
    
    out <- gsub("-", ":", gsub(",", "", gsub("chr", "", x)))
    out
}

#' list out genes
#' @param cnr a cnr bundle
#' 
#' @param at character, coordinate string e.g. "1:123456780:124567890", 
#' The \code{at} coordinates string is run through \code{convert_coord}, which
#' removes the preceeding `chr` and commas, and substitutes any dash `-`  with  `:` to
#' fit the format used here
#'
#'
#' @param identifier character, name of the gene identifier column to output,
#'  default "hgnc.symbol"
#' 
#' @return
#' A vector containing a list of genes or other gene identifier.
#'
#' @examples
#' data(cnr)
#' 
#' list_genes(cnr, at = "12:58000000:59000000")
#' 
#' @export
list_genes <- function(cnr, at, identifier = "hgnc.symbol") {
    cc <- unlist(strsplit(convert_coord(at), split = ":"))
    
    out1 <- cnr$gene.index[cnr$gene.index$chrom == cc[1], c(identifier, "chrom", "start", "end")]
    out1 <- out1[out1$start >= as.numeric(cc[2]), c(identifier, "chrom", "start", "end")]
    out1 <- out1[out1$end <= as.numeric(cc[3]), c(identifier, "chrom", "start", "end")]

    out <- out1[, identifier]

    return(out)
    
}
