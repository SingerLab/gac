#' merge lineage tracing results from Corey's Maximum Entropy method
#' @param cnr a cnr bundle
#' @param branch_table table named Branch_Table.csv
#' @param lineage_table table named Output_Table.csv
#' @param adjacency_alterations table named Tree_Edges.csv
#'
#' @examples \dontrun{
#' 
#' brc <- read.csv("data-raw/lineages/WD9048/Branch_Table.csv")
#' lin <- read.csv("data-raw/lineages/WD9048/Output_Table.csv")
#' adj <- read.csv("data-raw/lineages/WD9048/Tree_Edges.csv")
#' 
#' cnr <- merge_lineage_tables(WD9048a, brc, lin, adj)
#' 
#' }
#' 
#' @importFrom assertthat assert_that
#' 
#' @keywords internal
#' @noRd
merge_lineage_tables <- function(cnr, branch_table, lineage_table,
                                 adjacency_alterations) {

    assertthat::assert_that(all(branch_table$cellID %in% cnr$Y$cellID))
    assertthat::assert_that(nrow(branch_table) == nrow(cnr$Y))
    assertthat::assert_that(nrow(cnr$chromInfo) == nrow(lineage_table))
    assertthat::assert_that(all(rownames(cnr$chromInfo)== lineage_table$bin.id))
    
    cnr <- addPheno(cnr, df = branch_table)
    cnr <- addInfo(cnr, df = lineage_table)
    cnr$gene.index <- merge(cnr$gene.index, lineage_table, by = "bin.id")
    rownames(cnr$gene.index) <- cnr$gene.index$hgnc.symbol
    cnr[["adjacency_matrix"]] <- adjacency_alterations

    return(cnr)
} # end merge_lineage_tables
