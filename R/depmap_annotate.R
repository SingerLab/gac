#' annotate essential genes w/DepMap
#'
#' Takes a vector of depmap genes, and annotates the genes as esseintail
#' or not essential given a cutoff.  For details and downloads please visit
#' https://depmap.org/portal/ 
#'
#' @param cnr a cnr bundle
#' 
#' @param depMeans an object with the dependency score means for the cells
#' of interest
#'
#' @param cutoff cutoff of D2 score to consider essential gene or
#' non-essentail gene
#'
#' @return
#'
#' Function returns the gene.index with depmap annotation colums;
#'  both depMeans, and essential (essentian/non-essential)
#'
#' @examples
#'
#' \dontrun{
#' d2_map <- read.csv("inst/extdata/D2_combined_gene_dep_scores.csv.gz",
#'  header = TRUE, row.names = 1)
#'
#' ## select cells to use
#' keep <- grep("BONE", names(depmap))
#' 
#' d2_bone <- DepScores(depmap = d2_map, keep = keep)
#' 
#' depMeans <- rowMeans(d2_bone[,keep])
#' 
#' depmap_annotate(cnr, depMeans = depMeans)
#'
#' }
#' @keywords internal
depmap_annotate <- function(cnr, depMeans, cutoff = -0.5) {
## set boundaries at -0.8 for essential genes
    ngi <- cnr[["gene.index"]]
    ngi$depMeans <- depMeans[ngi$hgnc.symbol]
    ngi$essential <- "non-essential"
    ngi$essential[ngi$depMeans < cutoff] <- "essential"
    
    cnr[["gene.index"]] <- ngi
    
    return(cnr)
    
} ## depmap annotate

