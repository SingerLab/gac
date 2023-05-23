#' finds genes within one chromosome interval
#' 
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
#' list_genes_in_region(cnr, at = "12:58,000,000-59,000,000")
#'
#' list_genes_in_region(cnr, at = "12:58,000,000-59,000,000",
#'            show.columns = c("hgnc.symbol", "gene_biotype"))
#' 
#' 
#' @export
list_genes_in_region <- function(cnr, at,
                       show.columns = NULL,
                       identifier = "hgnc.symbol") {

    gene.list <- list_gene_symbols(cnr = cnr, at = at, identifier = identifier)
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
#' list_gene_symbols(cnr, at = "12:58000000:59000000")
#' 
#' @export
list_gene_symbols <- function(cnr, at, identifier = "hgnc.symbol") {
    cc <- unlist(strsplit(convert_coord(at), split = ":"))
    
    out1 <- cnr$gene.index[cnr$gene.index$chrom == cc[1], c(identifier, "chrom", "start", "end")]
    out1 <- out1[out1$start >= as.numeric(cc[2]), c(identifier, "chrom", "start", "end")]
    out1 <- out1[out1$end <= as.numeric(cc[3]), c(identifier, "chrom", "start", "end")]

    out <- out1[, identifier]

    return(out)
    
}


#' Pull gene details for a genomic region
#'
#' This function subsets the gene index for a genomic region of interest.
#'
#' @param cnr a cnr bundle
#'
#' @param chrom a chromosome name, must match one of 'cnr$gene.index$chrom'
#'
#' @param start region start
#'
#' @param end region end
#'
#' @return
#'
#' Returns the subset of the `gene.index` table for the genomic region.
#' 
#' @examples
#'
#' data(cnr)
#'
#' coord.df <- data.frame(chr = 12,
#'                     start = 69200804,
#'                     end = 69246466)
#'
#' get_gene_details(cnr, chrom = coord.df$chr,
#'                   start = coord.df$start, end = coord.df$end)
#'
#' coords.df <- data.frame(chr = c(1, 12), 
#'                      start = c(170120554, 69200804),
#'                      end =  c(172941951, 69246466))
#' do.call(rbind, apply(coords.df, 1, function(rr)
#'                      get_gene_details(cnr,
#'                                        chrom = rr[1],
#'                                        start = rr[2],
#'                                        end = rr[3])))
#' @keywords internal
#' @noRd
get_gene_details <- function(cnr, chrom = 12, start = 69200804, end = 69246466) {

    assertthat::assert_that(start < end)
    
    gene.details <- cnr$gene.index[as.character(cnr$gene.index$chrom) == chrom & cnr$gene.index$start > start & cnr$gene.index$end < end, ]

    return(gene.details)
    
} ## get_gene_details

    
#' Pull gene details for a genomic region
#'
#' This function subsets the gene index for a genomic region of interest.
#'
#' @param cnr a cnr bundle
#'
#' @param coord genomic region in ensembl format
#'
#' @return
#'
#' Returns the subset of the `gene.index` table for the genomic region.
#'
#' @examples
#'
#' data(cnr)
#' coord <- "12:69200804:69246466"
#' 
#' pull_gene_details(cnr, coord = "4:82351690:138565783")
#'
#' pull_gene_details(cnr, coord = coord)
#'
#' coords <- c("1:170120554:172941951",
#'           "12:69200804:69246466")
#' 
#' do.call(rbind, lapply(coords, function(rr)
#'                      pull_gene_details(cnr, coord = rr)))
#' 
#' @export
pull_gene_details <- function(cnr, coord = "12:69200804:69246466") {

    seqname <- unlist(strsplit(coord, split = ":"))[1]
    start <- as.numeric(unlist(strsplit(coord, split = ":"))[2])
    end <- as.numeric(unlist(strsplit(coord, split = ":"))[3])

    assertthat::assert_that(start < end)

    gene.details <- get_gene_details(cnr, chrom = seqname,
                                     start = start,
                                     end = end)

    return(gene.details)
} ## end pull gene details

