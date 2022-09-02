#' cluster heterogeneity across the sample
#'
#' @param cnr a cnr bundle
#'
#' @param by a categorical variable used to stratify the cell population.  If
#' NULL i.e. no stratification, the representation will be done overall
#'
#' @param cluster_column column to use as cluster
#'
#' @return
#' A `data.frame` with a simple cluster heterogeneity summary.  Table provides total
#'  counts of each clone (optionally be stratified by a group), clone representation,
#'  overal frequency per clone (also computed for stratified data), and Shannon's
#'  Diversity index for stratified data. 
#' 
#' @examples
#'
#' data(cnr)
#'
#' cnr <- phyloCNR(cnr)
#'
#' cnr <- setBrayClusters(cnr, tree.height = 0.065)
#'
#' cnr <- cluster_heterogeneity(cnr, cluster_column = "BrayC")
#'
#' @importFrom assertthat assert_that is.string
#' 
#' @export
cluster_heterogeneity <- function(cnr, by = NULL, cluster_column = NULL) {

    assertthat::assert_that(!is.null(cluster_column))
    assertthat::assert_that(any(cluster_column %in% names(cnr$Y)))
    
    if(is.null(by)) {

        occ <- table(cnr$Y[, cluster_column])
        
    } else {
        
        assertthat::is.string(by)
        assertthat::assert_that(by %in% names(cnr$Y))
        
        occ <- table(cnr$Y[, c(cluster_column, by)])
        
    }
    
    occ <- cluster_representation(occ)

    cnr$Y$final_cluster <- cnr$Y[, cluster_column]

    cnr[["cluster_heterogeneity"]] <- occ
    
    return(cnr)
} # cluster_heterogeneity


#' binary matrix 
#'
#' @param mat a matrix
#'
#' @export
binary.mat <- function(mat) {
    mat[mat >= 1] <- 1
    mat
} # end binary.mat


#' cluster representation
#'
#' @param cc cluster counts
#'
#' @examples
#' 
#' ## cluster table, no subsamples
#' cc <- table(sample(paste0("C", 1:12), size = 1000, replace = TRUE))
#' 
#' cc <- cluster_representation(cc)
#'
#' ## cluster table with multiple samples
#' cc <- data.frame(cluster = sample(paste0("C", 1:12), size = 1000,
#'                                   replace = TRUE),
#'             sample = sample(paste0("S", 1:6), size = 1000, replace = TRUE))
#' cc <- table(cc[, c("cluster", "sample")])
#' cc[cc <= 10] <- 0
#'
#' cc <- cluster_representation(cc)
#'
#' @importFrom vegan diversity
#' 
#' @keywords internal
#' @export
cluster_representation <- function(cc) {

    if(is.matrix(cc)) {
        n.cells.clone <- rowSums(cc)
        n.regions <- rowSums(binary.mat(cc))
        freqs <- cc / n.cells.clone
        overall.fq <- n.cells.clone / sum(n.cells.clone)
        colnames(freqs) <- paste0(colnames(cc), ".fq")
        sH <- apply(cc, 1, vegan::diversity) ## deault - shannon
        spatial.extent <- ifelse(n.regions >=2, "diffused", "local")
        
    } else {
        n.cells.clone <- sum(cc)
        n.regions <- 1
        overall.fq <- cc / n.cells.clone
        freqs <- cc / n.cells.clone
        sH <- NA
        spatial.extent <- NA
    }
    
    ccx <- data.frame(
        cbind(
            n.cells = n.cells.clone,
            overall.fq,
            n.regions = n.regions,
            sH = sH,
            spatial.extent = spatial.extent,
            cc,
            freqs))
    
    return(ccx)
    
} # end heterogeneity

