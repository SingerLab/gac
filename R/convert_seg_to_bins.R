#' convert .seg files to match a chromInfo matrix
#'
#' @param cnr a cnr bundle
#'
#' @param seg a table in .seg format 
#'
#' @param sample.id column name of sample ID in the .seg dat, default is "ID"
#'
#' @param bin.id column name containing a bin.id, defalt NULL,
#'  bin.id will be taken as the rownames of the data
#' 
#' @param coordinates.a  names of the columns containing chromosome, start,
#' and end coordinates in the cnr
#' 
#' @param coordinates.b  names of the columns containing chromosome, start,
#' and end coordinates in the seg data
#'
#' @param seg.mean.b name of column containg segment means default "seg.mean"
#'
#' @param multipcf weather to create a multi sample matrix of the .seg data
#'  or keep as run length encoding, akin to copynumber::multipcf,
#'  and copynumber::pcf, respectively.  This function requires
#'  tidyr::pivot_wider
#'
#' @examples
#' \dontrun{
#' data(cnr)
#' seg.data <- read.delim("path/to/data.seg")
#'
#' binned.seg <- seg2bins(cnr, seg.data,
#'                        sample.id = "sampleID",
#'                        bin.id = "bin.id")
#' }
#' 
#' @importFrom tidyr pivot_wider
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' 
#' @export
seg2bins <- function(cnr, seg, sample.id = "ID", bin.id = NULL,
                     coordinates.a = c("bin.chrom", "bin.start", "bin.end"),
                     coordinates.b = c("chrom", "loc.start", "loc.end"),
                     seg.mean.b = "seg.mean", multipcf = TRUE) {

    ## convert chromosomes to all numeric in both a and b
    numeric.chrom.a <- gsub("Y|y", "24",
                            gsub("X|x", "23", cnr$chromInfo[, coordinates.a[1]]))
    numeric.chrom.b <- gsub("Y|y", "24",
                            gsub("X|x", "23", seg[, coordinates.b[1]]))

    ## create genomic ranges objects for a and b
    a <- GenomicRanges::GRanges(seqnames = numeric.chrom.a,
           ranges = IRanges::IRanges(start = cnr$chromInfo[, coordinates.a[2]],
                                     end =  cnr$chromInfo[, coordinates.a[3]]))

    b <- GenomicRanges::GRanges(seqnames = numeric.chrom.b,
          ranges = IRanges::IRanges(start = seg[,coordinates.b[2]],
                                    end = seg[,coordinates.b[3]]))

    ## find overlaps between a and b
    overlap <- GenomicRanges::findOverlaps(b, a)

    ## if bin.id is given, use that column, otherwise it's calculated as
    ## 1:nrow of the cnr$chromInfo
    if(!is.null(bin.id)) {
        assertthat::assert_that(bin.id %in% colnames(cnr$chromInfo))

        ## expand segment data to bin information
        expanded.seg <- cbind(
            sampleID = seg[overlap@from,sample.id],
            bin.id = cnr$chromInfo[overlap@to, bin.id],
            cnr$chromInfo[overlap@to, coordinates.a],
            seg.mean = seg[overlap@from, seg.mean.b])
        
    } else {
        ## estimate bins from chromInfo
        bb <- 1:nrow(cnr$chromInfo)
        ## expand segment data to bin information
        expanded.seg <- cbind(
            sampleID = seg[overlap@from, sample.id],
            bin.id = bb[overlap@to],
            cnr$chromInfo[overlap@to, coordinates.a],
            seg.mean = seg[overlap@from, seg.mean.b])
        
    }
    ## multipcf: term imported from R/copynumber package
    ## pcf / multipcf are the functions to segment the data
    ## in multipcf, segments are broken down to be comon across
    ## all samples, similar to what a bin represents,
    ## multipcf format is in wide format,
    ## pcf is single-sample in long format
    if(multipcf) {
        ## pivot wider ; using median when a bin crosses over two segments
        pcf <- tidyr::pivot_wider(expanded.seg,
                                  id_cols = c("bin.id", all_of(coordinates.a)),
                                  names_from = "sampleID",
                                  values_from = "seg.mean",
                                  values_fn = median)
        ## merging with chromInfo to obtain all chromInfo bins
        ## removes bin.id
        pcf <- merge(cnr$chromInfo[, c(coordinates.a)],
                     pcf, all.x = TRUE)
        pcf <- pcf[, -c(grep("bin.id", names(pcf)))]
            
        
    } else {
        
        pcf <- expanded.seg
        
    }
    
    return(pcf)
    
} ## end seg2bins


#' run length encoding to create segment files from DDRC.df
#'
#' Still in development:  currently performs RLE by chromosome and
#' returns chromsome specific encoding.  Locantions are not currently 
#' offsest by chromsome start locations
#' 
#' @param cnr a cnr bundle
#'
#' @param chrom name of column in chromInfo to use as chromosomes
#'
#' @param pos name of column in chromInfo to use as chromosomal position
#'
#' @return
#' A run length encoding table from the clone profiles contained
#' in DDRC.df
#'
#'
#' @examples
#'
#' data(cnr)
#'
#' noisy.cells <- cnr$qc$cellID[cnr$qc$qc.status == "FAIL"]
#'
#' ## reduced pipeline to genrate DDRC clone profiles
#' cnr <- excludeCells(cnr, excl = noisy.cells)
#' cnr <- phylo_cnr(cnr, root.cell = "cell0")
#' cnr <- setBrayClusters(cnr)
#' cnr <- run_consensus_clustering(cnr, iters = 20, maxK = 40)
#' cnr <- doKSpectral(cnr)
#' cnr <- setKcc(cnr)
#' cnr <- cluster_heterogeneity(cnr, by = "category1",
#'           cluster_column = "ConsensusC")
#' cnr <- get_cluster_profiles(cnr)
#' 
#' ddrc_rle <- rle_ddrc(cnr)
#'
#' @importFrom assertthat assert_that
#' @importFrom reshape2 melt
#' 
#' @export
rle_ddrc <- function(cnr, chrom = "bin.chrom", pos = "bin.start") {
    
    assertthat::assert_that(!is.null(cnr[["DDRC.df"]]))

    resg <- cbind(cnr$chromInfo[, c(chrom, pos)],
                  cnr[["DDRC.df"]])
    resgm <- reshape2::melt(resg, id.vars = c(chrom, pos),
                            variable.name = "ID", value.name = "cn")
    resgm <- resgm[order(resgm$ID, resgm[,chrom]),]
    
    resgm$ID.CHROM <- paste(resgm$ID, resgm[,chrom], sep = "_")

    spl <- split(resgm$cn, f = resgm$ID.CHROM)

    nn <- names(spl)
    
    spl <- lapply(spl, rle)

    spl <- lapply(names(spl), function(ii) {
        rr = cbind(mean = spl[[ii]]$values,
                   nmark = spl[[ii]]$lengths)
        return(data.frame(rr))
    })
    
    for(i in 1:length(nn)) {
        spl[[i]] <- cbind(ID = gsub("(.*)_(.*)", "\\1", nn[[i]]),
                          chrom = gsub("(.*)_(.*)", "\\2", nn[[i]]),
                          spl[[i]],
                          start.bin = c(1, head(cumsum(spl[[i]]$nmark), -1)))
    }
    
    spl <- do.call(rbind, spl)
    spl <- spl[order(spl$ID, spl$chrom), ]
    spl$end.bin <- spl$start.bin + spl$nmark
    spl$end.bin[spl$start.bin == 1] <- spl$end.bin[spl$start.bin == 1] -1
    
    return(spl)
            
}
