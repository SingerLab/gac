#' proportion of polymorphic alleles as subclonal
#'
#' Estimate genome-wide proportion of polymorphic loci in a population
#'
#' @param cnr a cnr bundle
#'
#' @param exclude.chr chromosomes to exclude, default X, Y
#'
#' @param monomorphic.threshold alteration frequency at which a
#'  locus is to be considered monomorphic (universally altered).
#'  default 0.95
#'
#' @param noise.threshold lower-bound alteration frequency at
#'  which alterations are considered technical noise, rather
#'  than a true alteration, default 0.1
#'
#' @param chrom.col name of column containing chromosomes
#'
#' @return
#'
#' Single-value with the proportion of polymorphic loci for a given
#' sample population
#'
#' @examples
#' data(cnr)
#'
#' cnr <- get_alteration_frequencies(cnr)
#'
#' proportion_of_polymorphic_loci(cnr)
#'
#' @importFrom assertthat assert_that
#' @export
proportion_of_polymorphic_loci <- function(cnr, exclude.chr = c("X", "Y"),
                                           monomorphic.threshold = 0.95,
                                           noise.threshold = 0.1,
                                           chrom.col = "bin.chrom") {

    assertthat::assert_that(monomorphic.threshold > noise.threshold)
    

    ge <- subset_ci(cnr, exclude.chr = exclude.chr, chrom.col = chrom.col)
    
    n_pj <- sum(ge$altFQ <= monomorphic.threshold &
                ge$altFQ >= noise.threshold)
            
    n_total <- nrow(ge)

    P <- n_pj / n_total

    return(P)
}

#' Estimate the estimated number of alleles (copies) per locus
#'
#' Function is `experimental` and requires validation to asses it's value
#' 
#' @param cnr a cnr bundle
#'
#' @param exclude.chr chromosomes to exclude, default X, Y
#'
#' @param chrom.col name of column containing chromosomes
#'
#' @return
#' Estimates the number of alleles (copies per locus), across the population.
#'
#' Function is `experimental` and requires validation to asses it's value
#'
#' @examples
#' data(cnr)
#'
#' avg_num_alleles_per_locus(cnr)
#'
#' avg_num_alleles_per_locus(cnr, exclude.chr = c("X", "Y"))
#'
#' @importFrom assertthat assert_that
#' @export
avg_num_alleles_per_locus <- function(cnr, exclude.chr = NULL, chrom.col = "bin.chrom") {

    if(is.null(exclude.chr)) {

        K  <- nrow(cnr$X)
        ## list of allels per bin
        ni <- apply(cnr$X, 1, function(i) length(unique(as.numeric(i))))
        
    } else {
        assertthat::assert_that(all(exclude.chr %in% cnr$chromInfo[,chrom.col]))
    
        incl <- !which(cnr$chromInfo$bin.chrom %in% exclude.chr)
        K <- sum(incl)
        ## list of allels per bin
        ni <- apply(cnr$X[incl,], 1, function(i) length(unique(as.numeric(i))))
        
    }
        
    ## Average number of alleles per locus
    n = (1/K) * sum(ni)
    
    return(n)
}

#' calculate the proportion of the genome loss per cell
#' 
#' @param cnr a cnr bundle
#' 
#' @param by character, estimate percent genome loss at population, clone,
#'  or cell levels. Options are NULL, "clone", and "cell".  When NULL the
#'  estimate is done at the population level using bin alteration frequencies.
#'  Using "clone", and "cell", genome loss is performed using the clone profiles
#'  (DDRC.df), and cell profiles (X), respectively. Default NULL
#'
#' @param loss.threshold number of copies at which a locus is considered as a loss.
#'  Ignored when by = NULL. Default 1
#'
#' @param noise.threshold lower-bound alteration frequency at
#'  which alterations are considered technical noise, rather
#'  than a true alteration, default 0.1, ignored when by is "clone" or "cell"
#' 
#' @param genome.size size of the genome in base pairs.
#'  default 3098825702, human genome
#'
#' @param exclude.chr vector, chromosomes to exclude
#' 
#' @param chrom.col name of column containing chromosomes
#'
#' @return
#' The value of the percent genome loss across a population of cells.
#'
#' In the analysis by.clone = TRUE, the default is to consider all
#' copy numbers <= 1 as altered. For this, it's considered, at the moment,
#' that chromosomes X and Y have are to have two copies
#'
#' @examples
#'
#' data(cnr)
#'
#' noisy.cells <- cnr$qc$cellID[cnr$qc$qc.status == "FAIL"]
#'
#' ## reduced pipeline to genrate DDRC clone profiles
#' cnr <- excludeCells(cnr, excl = noisy.cells)
#' cnr <- phyloCNR(cnr, root.cell = "cell0")
#' cnr <- setBrayClusters(cnr)
#' cnr <- run_consensus_clustering(cnr, iters = 20, maxK = 40)
#' cnr <- doKSpectral(cnr)
#' cnr <- setKcc(cnr)
#' cnr <- cluster_heterogeneity(cnr, by = "category1",
#'           cluster_column = "ConsensusC")
#' cnr <- get_cluster_profiles(cnr)
#' cnr <- get_alteration_frequencies(cnr)
#'
#' percent_genome_loss(cnr)
#' 
#' percent_genome_loss(cnr, by = "cell")
#' 
#' percent_genome_loss(cnr, by = "clone")
#' 
#' @importFrom assertthat assert_that
#' 
#' @export
percent_genome_loss <- function(cnr,
                                by = NULL,
                                loss.threshold = 1,
                                noise.threshold = 0.1,
                                genome.size = 3098825702,
                                exclude.chr = NULL,
                                chrom.col = "bin.chrom") {

    ge <- subset_ci(cnr, exclude.chr = exclude.chr, chrom.col = chrom.col)
    
    if(is.null(by)) {
        assertthat::assert_that("delFQ" %in% names(cnr$chromInfo),
                                msg = "alteration frequencies not available in chromInfo.  Please run `get_alteration_frequencies`")
        
        pct.loss <- sum(ge[ge$delFQ >= noise.threshold, "bin.length"])  /
            genome.size
        
    } else {
        
        if(by == "clone") {

            assertthat::assert_that(!is.null(cnr$DDRC.df),
                                    msg = "`cnr$DDRC.df` not found")
            
            ddrc.loss <- cnr$DDRC.df <= loss.threshold
            
            pct.loss <- apply(ddrc.loss, 2, function(i) {
                sum(ge[i, "bin.length"]) / genome.size
            })
            
        } else {

            if(by == "cell") {

                cell.loss <- cnr$X <= loss.threshold
                pct.loss <- apply(cell.loss, 2, function(i) {
                    sum(ge[i, "bin.length"]) / genome.size
                })
            }
        }
    }
    
    return(pct.loss)
} ## end percent_genome_loss

#' calculate the proportion of the genome with gains
#' 
#' @param cnr a cnr bundle
#' 
#' @param by character, estimate percent genome loss at population, clone,
#'  or cell levels. Options are NULL, "clone", and "cell".  When NULL the
#'  estimate is done at the population level using bin alteration frequencies.
#'  Using "clone", and "cell", genome loss is performed using the clone profiles
#'  (DDRC.df), and cell profiles (X), respectively. Default NULL
#'
#' @param gain.thresholds lower and upper bound copy number at which a locus
#'  should considered a gain. The threshold is ignored when by = NULL.
#'  If only one value is provided, all copy numbers above that value are
#'  treated as a gain.  Default 3 to 5.
#'
#' @param noise.threshold lower-bound alteration frequency at
#'  which alterations are considered technical noise, rather
#'  than a true alteration, default 0.1, ignored when by is "clone" or "cell"
#' 
#' @param genome.size size of the genome in base pairs.
#'  default 3098825702, human genome
#'
#' @param exclude.chr vector, chromosomes to exclude
#' 
#' @param chrom.col name of column containing chromosomes
#' 
#' @return
#' The value of the percent genome gained across a population of cells.
#' gaines are considered when cell or clone copy numers are CN >= 3 & CN <= 5.
#' CN >= 6, it's considerd an amplification, consistent with \code{\link{callCN}}
#'
#' In the analysis by.clone = TRUE, the default is to consider all
#' copy numbers <= 1 as altered. For this, it's considered, at the moment,
#' that chromosomes X and Y have are to have two copies
#'
#' @examples
#'
#' data(cnr)
#'
#' noisy.cells <- cnr$qc$cellID[cnr$qc$qc.status == "FAIL"]
#'
#' ## reduced pipeline to genrate DDRC clone profiles
#' cnr <- excludeCells(cnr, excl = noisy.cells)
#' cnr <- phyloCNR(cnr, root.cell = "cell0")
#' cnr <- setBrayClusters(cnr)
#' cnr <- run_consensus_clustering(cnr, iters = 20, maxK = 40)
#' cnr <- doKSpectral(cnr)
#' cnr <- setKcc(cnr)
#' cnr <- cluster_heterogeneity(cnr, by = "category1",
#'           cluster_column = "ConsensusC")
#' cnr <- get_cluster_profiles(cnr)
#' cnr <- get_alteration_frequencies(cnr)
#'
#' percent_genome_gain(cnr)
#' 
#' percent_genome_gain(cnr, by = "cell")
#' 
#' percent_genome_gain(cnr, by = "clone")
#' 
#' @importFrom assertthat assert_that
#' 
#' @export
percent_genome_gain <- function(cnr,
                                by = NULL,
                                gain.thresholds = c(3, 5),
                                noise.threshold = 0.1,
                                genome.size = 3098825702,
                                exclude.chr = NULL,
                                chrom.col = "bin.chrom") {

    
    ge <- subset_ci(cnr, exclude.chr = exclude.chr, chrom.col = chrom.col)
    
    if(is.null(by)) {

        assertthat::assert_that("AmpFQ" %in% names(cnr$chromInfo),
                                msg = "alteration frequencies not available in chromInfo.  Please run `get_alteration_frequencies`")
        
        pct.gain <- sum(ge[ge$AmpFQ >= noise.threshold,"bin.length"]) / genome.size
        
    } else {
        
        if(by == "clone") {
            
            assertthat::assert_that(!is.null(cnr$DDRC.df),
                                    msg = "`cnr$DDRC.df` not found")
            
            if(length(gain.thresholds) == 1) {
                ddrc.gain <- cnr$DDRC.df >= gain.thresholds

            } else {

                ddrc.gain <- cnr$DDRC.df >= gain.thresholds[1] &
                    cnr$DDRC.df <= gain.thresholds[2]

            }

            pct.gain <- apply(ddrc.gain, 2, function(i) {
                sum(ge[i, "bin.length"]) / genome.size

            })
            
        } else {

            if(by == "cell")

                if(length(gain.thresholds) == 1) {

                    cell.gain <- cnr$X >= gain.thresholds

                } else {

                    cell.gain <- cnr$X >= gain.thresholds[1] &
                        cnr$X <= gain.thresholds[2]

                }    

            pct.gain <- apply(cell.gain, 2, function(i) {
                sum(ge[i, "bin.length"]) / genome.size

            })
        }
    }

    return(pct.gain)
    
} ## end percent_genome_gain

#' calculate the proportion of the genome amplified
#' @param cnr a cnr bundle
#' 
#' @param by character, estimate percent genome loss at population, clone,
#'  or cell levels. Options are "clone", and "cell".
#'  Using "clone", and "cell", genome loss is performed using the clone profiles
#'  (DDRC.df), and cell profiles (X), respectively. Default "cell"
#'
#' @param amplification.threshold number of copies at which a locus is considered as an amplification.
#'
#' @param ... additional parameters to \code{\link{percent_genome_gain}}
#'
#' @examples
#'
#' data(cnr)
#'
#' noisy.cells <- cnr$qc$cellID[cnr$qc$qc.status == "FAIL"]
#'
#' ## reduced pipeline to genrate DDRC clone profiles
#' cnr <- excludeCells(cnr, excl = noisy.cells)
#' cnr <- phyloCNR(cnr, root.cell = "cell0")
#' cnr <- setBrayClusters(cnr)
#' cnr <- run_consensus_clustering(cnr, iters = 20, maxK = 40)
#' cnr <- doKSpectral(cnr)
#' cnr <- setKcc(cnr)
#' cnr <- cluster_heterogeneity(cnr, by = "category1",
#'           cluster_column = "ConsensusC")
#' cnr <- get_cluster_profiles(cnr)
#' cnr <- get_alteration_frequencies(cnr)
#' 
#' percent_genome_amplified(cnr)
#' 
#' percent_genome_amplified(cnr, by = "cell")
#' 
#' percent_genome_amplified(cnr, by = "clone")
#' 
#' @export
percent_genome_amplified <- function(cnr, by = "cell",
                                     amplification.threshold = 6,
                                     ...) {
    if(is.null(by)) {
        warning('output is equivalent as running percent_genome_gain
                  with by = NULL')
        pct.amp <- percent_genome_gain(cnr, by = by, ...)
    } else {
        pct.amp <- percent_genome_gain(
            cnr, by = by,
            gain.thresholds = amplification.threshold,
            ...)
    }
    return(pct.amp)
} ## end percent_genome_amplified


#' subset for chromosomes of interest
#' @param cnr a cnr bundle
#' 
#' @param exclude.chr vector, chromosomes to exclude
#' 
#' @param chrom.col name of column containing chromosomes
#'
#' @return
#' conditional subset of chromInfo w/o excluded chromosomes
#' 
#' @keywords internal
subset_ci <- function(cnr, exclude.chr = NULL, chrom.col = "bin.chrom") {

    if(is.null(exclude.chr)) {
        
        ge <- cnr$chromInfo
        
    } else {
        
        assertthat::assert_that(chrom.col %in% names(cnr$chromInfo))
        assertthat::assert_that(all(exclude.chr %in% cnr$chromInfo[,chrom.col]))
        ge <- cnr$chromInfo[!cnr$chromInfo[,chrom.col] %in% exclude.chr, ]

    }
    
    return(ge)
    
}
