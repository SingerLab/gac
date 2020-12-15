#' Pull gene details for a genomic region
#'
#' This function subsets the gene index for a genomic region of interest.
#'
#' @param cnr a cnr bundle
#'
#' @param seqnames a chromosome name, must match 'cnr$gene.index$seqnames'
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
#' library(dplyr)
#' data(cnr)
#' coord.df <- data.frame(chr = 12,
#'                     start = 69200804,
#'                     end = 69246466)
#'
#' get_gene_details(cnr, seqnames = coord.df$chr,
#'                   start = coord.df$start, end = coord.df$end)
#'
#' coords.df <- data.frame(chr = c(1, 12), 
#'                      start = c(170120554, 69200804),
#'                      end =  c(172941951, 69246466))
#' do.call(rbind, apply(coords.df, 1, function(rr)
#'                      get_gene_details(cnr,
#'                                        seqnames = rr[1],
#'                                        start = rr[2],
#'                                        end = rr[3])))
#' 
#' @export
get_gene_details <- function(cnr, seqnames = 12, start = 69200804, end = 69246466) {

    assertthat::assert_that(start < end)
    
    gene.details <- cnr$gene.index[as.character(cnr$gene.index$seqnames) == seqnames & cnr$gene.index$start > start & cnr$gene.index$end < end, ]

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
#' library(dplyr)
#' data(cnr)
#' coord <- "12:69200804:69246466"
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
    start <- unlist(strsplit(coord, split = ":"))[2]
    end <- unlist(strsplit(coord, split = ":"))[3]

    assertthat::assert_that(start < end)

    gene.details <- get_gene_details(cnr, seqnames = seqname,
                                     start = start,
                                     end = end)

    return(gene.details)
} ## end pull gene details

