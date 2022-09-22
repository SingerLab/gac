#' Heatmap of cells with VDJ recombination
#'
#' @param cnr a cnir bundle
#'
#' @param vdj.genes list of VDJ genes.  By default is set to NULL, vdj.genes is
#'  the union of T-cell Receptor and B-cell Receptor genes (i.e. TR, IG genes,
#'  used in genotype_vdj)
#'                  
#'
#' @param vdjGeneAnno HeatmapAnnotation object to annotate gene attributes
#'
#' @param vdjCellAnno rowAnnotation object to annotate cells
#'
#' @param ... additional arguments passed to ComplexHeatmap::Heatmap
#'
#' @return
#' Returns a heatmap of VDJ cells with rows as cells, and TR/IG genes as columns.
#'  Heatmap is split by chromosome name, and cell types are annotated on the left
#' 
#' @examples
#' data(cnr)
#'
#' cnr <- genotype_vdj(cnr)
#'
#' vdjHeatmap(cnr)
#'
#' @importFrom ComplexHeatmap HeatmapAnnotation rowAnnotation Heatmap
#' @export
vdjHeatmap <- function(cnr, vdj.genes = NULL, vdjGeneAnno = NULL, 
                       vdjCellAnno = NULL, ...) {
    
    vdjCells <- names(cnr$vdj.cells)

    cnr$gene.index$chrom <- factor(cnr$gene.index$chrom, levels = c(1:22, "X", "Y"))
    
    if(is.null(vdj.genes)) {
        ## get T-cell receptor VDJ genes from gene.index
        tr.genes <- cnr$gene.index$hgnc.symbol[grepl("TR.*gene", cnr$gene.index$gene.type)]
        tr.genes <- tr.genes[tr.genes %in% colnames(cnr$genes)]
        tr.genes <- tr.genes[!grepl("pseudogene", cnr$gene.index[tr.genes, "gene.type"])]
        ## get B-cell receptor VDJ genes from gene.index
        ig.genes <-  cnr$gene.index$hgnc.symbol[grepl("IG.*gene", cnr$gene.index$gene.type)]
        ig.genes <- ig.genes[ig.genes %in% names(cnr$genes)]
        ig.genes <- ig.genes[!grepl("pseudogene", cnr$gene.index[ig.genes, "gene.type"])]
        ## vdj genes
        vdj.genes <- union(tr.genes, ig.genes)
    } else {
        vdj.genes <- vdj.genes
    }
    
    if(is.null(vdjGeneAnno)) {
        vdjChr <- ComplexHeatmap::HeatmapAnnotation(
            "Chr" = cnr$gene.index[gsub("\\.", "-", vdj.genes), "chrom"],
            "Gene Type" = cnr$gene.index[gsub("\\.", "-", vdj.genes), "gene.type"],
            annotation_name_side = "left",
            show_legend = TRUE)
    } else {
        vdjChr <- vdjGeneAnno
    }

    if(is.null(vdjCellAnno)) {
        vdjAnnot <- ComplexHeatmap::rowAnnotation(
            "Cell Type" = cnr$vdj.cells,
            show_legend = TRUE)
    } else {
        vdjAnnot <- vdjCellAnno
    }

    vdjH <- ComplexHeatmap::Heatmap(
                    cnr$genes[vdjCells, vdj.genes],
                    cluster_rows = TRUE,
                    cluster_columns = FALSE,
                    cluster_row_slices = FALSE,
                    top_annotation = vdjChr,
                    left_annotation = vdjAnnot,
                    column_split = cnr$gene.index[gsub("\\.", "-", vdj.genes), "chrom"],
                    ...)
    
    vdjH

}
