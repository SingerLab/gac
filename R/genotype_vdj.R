#' Genotype cells at the VDJ recombination sites
#'
#' @param cnr a cnr bundle.  gene.index MUST include a column named `gene.type`
#'   with pattern "TR.*gene" for T-cell receptor genes, and "IG.*gene" for B-cell
#'   receptor genes.  These are the defaults for gene.type in bioconductor::biomaRt
#'
#' @param gene.type.column name of gene type column
#'
#' @return
#' Returns a CNR object with `cell.type` column containing assigments for
#'  `T-cell`, `B-cell`, and `vdj.unspecified`.  The function will only annotate
#'  cells with deletions on a `TR` or `IG` gene present in the gene.index table.
#'  All other annotation will be left intact.  If `cell.type` column is not
#'  present one will be created with `NA`.
#'
#' An additional cnr[["vdj.cells"]] vector is also returned containing the cell
#'  type of only VDJ cells.
#'
#' @examples
#' data(cnr)
#'
#' cnr <- genotype_vdj(cnr)
#' 
#' @importFrom assertthat assert_that
#' 
#' @export
genotype_vdj <- function(cnr, gene.type.column = "gene.type") {

    assertthat::assert_that(gene.type.column %in% names(cnr$gene.index))

    if(!"cell.type" %in% names(cnr$Y)) {
        cnr$Y$cell.type <- "unspecified.cell"
    }

    if("cell.type" %in% names(cnr$Y)) {
        assertthat::assert_that(is.character(cnr$Y$cell.type))
    } else {
        cnr$Y$cell.type <- as.character(cnr$Y$cell.type)
    }
    
    ## get T-cell receptor VDJ genes from gene.index
    tr.genes <- cnr$gene.index$hgnc.symbol[grepl("TR.*gene", cnr$gene.index$gene.type)]
    tr.genes <- tr.genes[tr.genes %in% colnames(cnr$genes)]
    tr.genes <- tr.genes[!grepl("pseudogene", cnr$gene.index[tr.genes, gene.type.column])]

    ## get B-cell receptor VDJ genes from gene.index
    ig.genes <-  cnr$gene.index$hgnc.symbol[grepl("IG.*gene", cnr$gene.index$gene.type)]
    ig.genes <- ig.genes[ig.genes %in% names(cnr$genes)]
    ig.genes <- ig.genes[!grepl("pseudogene", cnr$gene.index[ig.genes, gene.type.column])]

    ## get copy number states for IG and TR genes
    trGeno <-  cnr$genes[, tr.genes]
    igGeno <-  cnr$genes[, ig.genes]

    ## look for cells with loss or deletion of TR genes
    trCells <- apply(trGeno, 1, function(i) sum(i >= 2) != length(tr.genes))
    igCells <- apply(igGeno, 1, function(i) sum(i >= 2) != length(ig.genes))

    ## get cell names
    vdjCells <- union(names(igCells)[igCells], names(trCells)[trCells])
    vdjType <- rep(NA, length(vdjCells))

    ## assing cell type:  if deletions occur at IG-genes, assing B-cell
    ## assing cell type:  if deletions occur at TR-genes, assign T-cell
    ## if on both, assing vdj.unspecified
    vdjType[vdjCells %in% names(igCells)[igCells]] <- "B-cell"
    vdjType[vdjCells %in% names(trCells)[trCells]] <- "T-cell"
    vdjType[vdjCells %in% intersect(names(igCells)[igCells], names(trCells)[trCells])] <- "vdj.unspecified"
    names(vdjType) <- vdjCells

    assertthat::assert_that(all(names(vdjCells) == names(cnr$X)))
    assertthat::assert_that(all(names(vdjCells) == rownames(cnr$Y)))
    
    cnr$Y[cnr$Y$cellID %in% names(vdjType)[vdjType == "T-cell"], "cell.type"] <- "T-cell"
    cnr$Y[cnr$Y$cellID %in% names(vdjType)[vdjType == "B-cell"], "cell.type"] <- "B-cell"
    cnr$Y[cnr$Y$cellID %in% names(vdjType)[vdjType == "vdj.unspecified"], "cell.type"] <- "vdj.unspecified"

    cnr[["vdj.cells"]] <- vdjType

    return(cnr)
    
} ## end genotype_vdj
