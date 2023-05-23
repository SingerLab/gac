#' Heatmap of cells with VDJ recombination
#'
#' @param cnr a cnir bundle
#'
#' @param vdj.genes list of VDJ genes.  By default is set to NULL, vdj.genes is
#'  the union of T-cell Receptor and B-cell Receptor genes (i.e. TR, IG genes,
#'  used in genotype_vdj)
#'                  
#' @param vdjGeneAnno HeatmapAnnotation object to annotate gene attributes
#'
#' @param vdjCellAnno rowAnnotation object to annotate cells
#'
#' @param gene.type.column name of column with gene type, default "gene_biotype
#'
#' @param col copy number colors, default is NULL and will generate a black, white
#'  and red color palette.
#'
#' @param cell.type.colors color map for cell types. default is NULL
#'
#' @param chromosome.colors color map for VDJ chromosomes. default is NULL
#'
#' @param gene.type.colors color map for VDJ Gene types. default is NULL
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
                       vdjCellAnno = NULL, gene.type.column = "gene_biotype",
                       col = NULL, cell.type.colors = NULL,
                       chromosome.colors =NULL, gene.type.colors = NULL,
                       ...) {

    if(is.null(col)) {
        vdj.colors <- circlize::colorRamp2(
                                    breaks = c(0, 1, 2, 3),
                                    c("#242424", "#242424", "#FFFFFF", "#FF523F"))
    } else {
        vdj.colors <- col
    }
    
    vdjCells <- names(cnr$vdj.cells)

    if(!is.factor(cnr$gene.index$chrom)) {
        cnr <- order_genes(cnr)
    }

    ## default list of vdj genes based on gene.index
    if(is.null(vdj.genes)) {
        ## get T-cell receptor VDJ genes from gene.index
        tr.genes <-
            cnr$gene.index$hgnc.symbol[grepl("TR.*gene",
                                             cnr$gene.index[, gene.type.column])]
        tr.genes <- tr.genes[tr.genes %in% colnames(cnr$genes)]
        tr.genes <- tr.genes[!grepl("pseudogene",
                                    cnr$gene.index[tr.genes, gene.type.column])]
        ## get B-cell receptor VDJ genes from gene.index
        ig.genes <-
            cnr$gene.index$hgnc.symbol[grepl("IG.*gene",
                                             cnr$gene.index[,gene.type.column])]
        ig.genes <- ig.genes[ig.genes %in% names(cnr$genes)]
        ig.genes <- ig.genes[!grepl("pseudogene",
                                    cnr$gene.index[ig.genes, gene.type.column])]
        ## vdj genes
        vdj.genes <- union(tr.genes, ig.genes)
    } else {
        vdj.genes <- vdj.genes
    }

    ## build VDJ gene annotation
    if(is.null(gene.type.colors)) {
        gene.type.colors <- c(IG_C_gene = "#7FC97F", IG_J_gene = "#BEAED4",
                              TR_C_gene = "#FDC086", TR_J_gene = "#FFFF99",
                              TR_V_gene = "#386CB0", TR_D_gene = "#F0027F")
    }
    
    if(is.null(chromosome.colors)) {
        chromosome.colors <- c("2" = "#404040", "7" = "#BABABA",
                               "14" = "#404040", "22" = "#BABABA")
    }
    

    if(is.null(vdjGeneAnno)) {
        vdjChr <- ComplexHeatmap::HeatmapAnnotation(
            "Chr" = cnr$gene.index[gsub("\\.", "-", vdj.genes), "chrom"],
            "Gene Type" = cnr$gene.index[gsub("\\.", "-", vdj.genes),
                                         gene.type.column],
            annotation_name_side = "left",
            col = list(
                Chr = chromosome.colors,
                "Gene Type" = gene.type.colors),
            show_legend = TRUE)
    } else {
        vdjChr <- vdjGeneAnno
    }

    ## build VDJ cell annotation
    if(is.null(cell.type.colors)) {
        cell.type.colors <- c(
            "B-cell" = "#F1A340",
            "T-cell" = "#33A02C", 
            "vdj.unspecified" = "#CAB2D6")
    }
    
    if(is.null(vdjCellAnno)) {
        vdjAnnot <- ComplexHeatmap::rowAnnotation(
            "Cell Type" = cnr$vdj.cells,
            col = list(
                "Cell Type" = cell.type.colors),
            show_legend = TRUE)
    } else {
        vdjAnnot <- vdjCellAnno
    }

    vdjH <- ComplexHeatmap::Heatmap(
                    cnr$genes[vdjCells, vdj.genes],
                    col = vdj.colors,
                    cluster_rows = TRUE,
                    cluster_columns = FALSE,
                    cluster_row_slices = FALSE,
                    top_annotation = vdjChr,
                    left_annotation = vdjAnnot,
                    column_split = cnr$gene.index[gsub("\\.", "-", vdj.genes),
                                                  "chrom"],
                    ...)
    
    vdjH

}
