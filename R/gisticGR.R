#' Gistic Genomic Ranges, get genes in gene.index that overlap with significant
#' genomic ranges in GISTIC2 output
#'
#' **Still in development**
#' ** requires GenommicRanges
#'
#' note if oncokb or depmeans are FALSE, both annotations will be ignored
#' 
#' @param cnr a cnr bundle
#'
#' @param grTR gistic regions table, output from gisticRegions
#'
#' @param oncokb OncoKB cancer gene annotations (TRUE/FALSE)
#'   (see ?oncoKB_annotation)
#'
#' @param depmeans dependency score means for cells of interest
#'    e.g. DepMap achilles scores (see ?depmap_annotation), (TRUE/FALSE)
#'
#' @importFrom assertthat assert_that
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges findOverlaps
#' 
#' @keywords internal
#' @noRd
gisticGR <- function(cnr, grTR, oncokb = TRUE, depmeans = TRUE) {
    ## set up vector of amplified regions
    amp <- grTR[grep("Amp", grTR$Unique.Name), c("chr", "start", "end")]
    ## set up vector of deleted regions
    del <- grTR[grep("Del", grTR$Unique.Name), c("chr", "start", "end")]
    
    gr <- GenomicRanges::GRanges(seqnames = c(amp$chr, del$chr),
                  ranges = IRanges::IRanges(start = as.numeric(c(amp$start, del$start)), end = as.numeric(c(amp$end, del$end))),
                  alteration.type = c(rep("Amplification", nrow(amp)), rep("Deletion", nrow(del))))

    if(oncokb & depmeans) {
        assertthat::assert_that("oncoKB" %in% names(cnr$gene.index))
        assertthat::assert_that("depMeans" %in% names(cnr$gene.index))
        assertthat::assert_that("essential" %in% names(cnr$gene.index))

        geneIndex <- GenomicRanges::GRanges(seqnames = cnr$gene.index$chrom,
                         ranges = IRanges::IRanges(start = cnr$gene.index$start,
                                          end = cnr$gene.index$end),
                         strand = cnr$gene.index$strand,
                         ensembl_gene_id = cnr$gene.index$ensembl_gene_id,
                         hgnc.symbol = cnr$gene.index$hgnc.symbol,
                         gene.type = cnr$gene.index$gene.type,
                         bind.id = cnr$gene.index$bin.id,
                         oncoKB = cnr$gene.index$oncoKB,
                         depMeans = cnr$gene.index$depMeans,
                         essential = cnr$gene.index$essential)
    
    overlaps <- GenomicRanges::findOverlaps(geneIndex, gr)
    
    gistic.genes <- data.frame(gr[overlaps@to,], geneIndex[overlaps@from,])
    
    names(gistic.genes) <- c("chrom", "start", "end", "width", "strand",
                             "alteration.type", 
                             "chr", "gene.start", "gene.end", "gene.width",
                             "gene.strand", "ensembl_gene_id",
                             "hgnc.symbol", "gene.type", "bin.id",
                             "oncoKB", "depMeans", "essential")

    } else {
        
        geneIndex <-
            GenomicRanges::GRanges(
                               seqnames = cnr$gene.index$chrom,
                               ranges = IRanges::IRanges(
                                                     start = cnr$gene.index$start,
                                                     end = cnr$gene.index$end),
                               strand = cnr$gene.index$strand,
                               ensembl_gene_id = cnr$gene.index$ensembl_gene_id,
                               hgnc.symbol = cnr$gene.index$hgnc.symbol,
                               gene.type = cnr$gene.index$gene.type,
                               bind.id = cnr$gene.index$bin.id)
        
        overlaps <- GenomicRanges::findOverlaps(geneIndex, gr)
        
        gistic.genes <- data.frame(gr[overlaps@to,], geneIndex[overlaps@from,])
        
        names(gistic.genes) <- c("chrom", "start", "end", "width", "strand",
                                 "alteration.type", 
                                 "chr", "gene.start", "gene.end", "gene.width",
                                 "gene.strand", "ensembl_gene_id",
                                 "hgnc.symbol", "gene.type", "bin.id")
        
    }

    return(gistic.genes)
    
} ## end gisticGR

