#' Process GISTIC2 output // output only grTR
#' read in gistic leasions file
#' provide path to gistic directory w/ all_lesions.conf_**.txt
#'
#' **Still in development**
#'
#' @param gisticDir path to directory of gistic output
#'
#' @param conf confidence used in gistic, default is 75, for single-cells
#' 80 or highrer works a little better.  Also, requires increasing amp/del
#' calling thresholds to log2(1/2) for deletions, and log2(3/2) for gains
#' and amplifications
#'
#' @examples
#' \dontrun{
#'
#' grTR <- gisticRegions("path/to/gistic/", conf = 80)
#'
#' }
#' 
#' @keywords internal
#' @noRd
gisticRegions <- function(gisticDir, conf = 80) {
    gr <- read.delim(file.path(gisticDir, paste0("all_lesions.conf_", conf, ".txt")),
                     stringsAsFactors = FALSE)
    
    ## parse out genomic coordinates, and bin coordinates
    gr$chr <- gsub("chr", "", gsub("(chr[0-9XY]+):.*", "\\1", gr$Wide.Peak.Limits))
    gr$start <- gsub("(chr[0-9XY]+):([0-9]+)-([0-9]+)\\(.*", "\\2", gr$Wide.Peak.Limits)
    gr$end <- gsub("(chr[0-9XY]+):([0-9]+)-([0-9]+)\\(.*", "\\3", gr$Wide.Peak.Limits)
    gr$bin.start <- gsub("(chr[0-9XY]+):([0-9]+)-([0-9]+)\\(probes ([0-9]+):([0-9]+)\\) *", "\\4", gr$Wide.Peak.Limits)
    gr$bin.end <- gsub("(chr[0-9XY]+):([0-9]+)-([0-9]+)\\(probes ([0-9]+):([0-9]+)\\) *", "\\5", gr$Wide.Peak.Limits)
    
    ## get threshold and copy number data
    grCN <- gr[grep("Actual Copy Change Given", gr$Amplitude.Threshold),]
    grTR <- gr[grep("Actual Copy Change Given", gr$Amplitude.Threshold,
                    invert = TRUE),]

    return(grTR)
    
} ## gisticRegions

