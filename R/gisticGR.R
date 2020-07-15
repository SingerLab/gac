#' Gistic Genomic Ranges, get genes in gene.index that overlap with significant genomic ranges in GISTIC2 output
#'
#' **Still in development**
#' ** requires GenommicRanges
#' 
#' @param cnr a cnr bundle
#'
#' @param grTR gistic regions table, output from gisticRegions
#'
#' @param oncokb OncoKB annotation (downloaded from OncoKB)
#'
#' @param depmeans dependency score means for cells of interest
#'    e.g. DepMap achilles scores
#'
#' @param cutoff threshold to consider a gene as essential or non-essential
#'
#' @import GenomicRanges
#' 
#' @export
gisticGR <- function(cnr, grTR, oncokb, depmeans, cutoff = -0.5) {
    ## set up vector of amplified regions
    amp <- grTR[grep("Amp", grTR$Unique.Name), c("chr", "start", "end")]
    ## set up vector of deleted regions
    del <- grTR[grep("Del", grTR$Unique.Name), c("chr", "start", "end")]
    
    gr <- GenomicRanges::GRanges(seqnames = c(amp$chr, del$chr),
                  ranges = IRanges::IRanges(start = as.numeric(c(amp$start, del$start)), end = as.numeric(c(amp$end, del$end))),
                  alteration.type = c(rep("Amplification", nrow(amp)), rep("Deletion", nrow(del))))
    
    geneIndex <- GenomicRanges::GRanges(seqnames = cnr$gene.index$seqnames,
                         ranges = IRanges::IRanges(start = cnr$gene.index$start,
                                          end = cnr$gene.index$end),
                         strand = cnr$gene.index$strand,
                         ensembl_gene_id = cnr$gene.index$ensembl_gene_id,
                         hgnc.symbol = cnr$gene.index$hgnc.symbol,
                         gene.type = cnr$gene.index$gene.type,
                         bind.id = cnr$gene.index$bin.id)
    
    overlaps <- GenomicRanges::findOverlaps(geneIndex, gr)
    
    gistic.genes <- data.frame(gr[overlaps@to,], geneIndex[overlaps@from,])
    
    names(gistic.genes) <- c("seqnames", "chrom", "start", "end", "strand",
                             "alteration.type", 
                             "chr", "gene.start", "gene.end", "width",
                             "gene.strand", "ensembl_gene_id",
                             "hgnc.symbol", "gene.type", "bin.id")

    gistic.genes$oncokb <- "not.oncokb"

    gistic.genes$oncokb[gistic.genes$hgnc.symbol %in% oncokb$Hugo.Symbol] <- "cancer.gene"
    gistic.genes$oncokb[gistic.genes$hgnc.symbol %in% oncokb$Hugo.Symbol[oncokb$Is.Oncogene]] <- "oncogene"
    gistic.genes$oncokb[gistic.genes$hgnc.symbol %in% oncokb$Hugo.Symbol[oncokb$Is.Tumor.Supressor.Gene]] <- "tsg"
    
    gistic.genes$depMeans <- depmeans[gistic.genes$hgnc.symbol]
    gistic.genes$essential <- "non-essential"
    gistic.genes$essential[gistic.genes$depMeans < cutoff] <- "essential"
    
    return(gistic.genes)
    
}

