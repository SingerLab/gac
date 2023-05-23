#'  cnr: copy number, rounded
#'
#' @description
#' 
#' The cnr object is a list of four relational matrices. The bins, genes,
#' annotation, qc, chromInfo, and gene.index.  The structure is inspired by
#' Scanpy's AnnData which cleverly integrates complex data into a simple
#' architecture.
#'
#' The object stores processing results required for data exploration,
#' visualization and genetic analyses.  Specifically, the distance matrix,
#' phylogenetic analysis, pseudobulk, and heterogeneity analysis.
#'
#'
#' @format An object class list containg a rounded CNR
#'
#' * Input 
#' \itemize{
#' 
#'   \item X, An integer matrix of bins x n.cells containing copy number
#' estimations for each bin[i] and cell[j] Where bins represent a common genomic
#' segment across all cells (either a fixed with, variable binning, or .seg data). 
#'  This data can be constructed using a variable length bin and CBS (Varbin
#' algorithm) (Baslan et al 2012.), and implemented on Ginkgo, or from hmmCopy.
#' These are upstream analyses to the package.
#'
#' For mutations in single-cells, the cnq can be a binary incidence (0,1) matrix
#'  representing presence or absence of specific mutations, or a ternary (0,1,2)
#'  representing genotypes as the number of alternate allele copies
#'
#' We noticed that estimates for deletions and losses don't follow standard
#'  numeric rounding. We set a the thresholds of < 0.2 (average quantal estimate
#'  of the Y chromosome in females) for deletions (i.e. 0 copies); between 0.2
#'  and 1.2 for losses (i.e. 1 copy); between 1.2 and 2.5 for 2 copies, and
#'  everything else standard numeric rounding.
#' 
#'   \item genes, gene copy number interpolation from bins.  The genes matrix
#' is an interpolated, transposed, expansion of bins. The expansion is
#' constructed internally using the expand2genes function.
#' 
#'   \item Y, phenotypic data of single-cells, contains cells as rows, and
#' phentypes in columns. Phenotypes can be information about individual samples,
#' or if same-cell methods were used, the RNA expression from the same cell.  
#'
#'   \item qc, quality control metrics. This matrix contains additional metadata
#' that is technical, e.g. number of reads, MAPD estimates, and the PASS/FAIL
#' qc.status for individual cells.  contains cells as rows and metadata as columns
#'
#'  \item chromInfo, ordered chromosome information for the bins
#'
#'  \item gene.index, table to map bins to genes
#'
#'  \item cells, a list of cells
#'
#'  \item bulk, logical, weather the data is bulk DNA or cells.  If TRUE, data
#' will not be rounded and it's assumed is log ratio data.  If FALSE, data is
#' considered as single-cell and copy numbers are considered as integer copy
#' number.  Estimates are rounded to the nearest integer (see above).
#'   ...
#' }
#'
#' * Output
#' \itemize{
#'   \item cdb, pairwise cell dissimilarity using Bray-Curtis
#'
#'   \item hcdb, heirarchical clustering of cells based
#' 
#'   \item phylo, cell phyogenetic tree.  Analysis is produced with
#'     \code{\link[ape]}. Default is "balanced minimum evolution"
#' 
#'   \item tree.height, height cutoff of the tree, set as intersection between
#'      the total number of multi-cell clusters and one-cell clusters
#' 
#'   \item ccp, output from ConsensusClusterPlus, the first element is a color map.
#'      Elements 2:(40) contain 4 outputs: consensusMatrix, consensusTree,
#'      consensusClass, ml, and clrs.  Each element corresponds to k in k-means
#'      clustering
#' 
#'   \item kStats, spectral analysis of consensus clustering
#' 
#'   \item eigenVals, spectral analysis eigen values
#' 
#'   \item optK optimumn k-parameter (kCC), and stable K (sK)
#' 
#'   \item cluster_heterogeneity, summary metrics of clones/clusters
#' 
#'   \item uclust, unique clusters with 3 or more cells
#' 
#'   \item DDRC.df, bin pseudobulk profiles of each clone (or other chosen group)
#' 
#'   \item DDRC.g, gene pseudobulk profiles
#' 
#'   \item vdj.cells list of vdj.cells produced by \code{genotype_vdj}
#' 
#'   \item DDRC.dist, bray curtis disimilarity of clones
#' 
#'   \item DDRC.phylo, phylogenetic analysis of clones
#' 
#' }
#' 
#' @docType data
#' @keywords datasets
#' @usage data(cnr)
#' @source \url{https://github.com/SingerLab/gac"}
"cnr"
