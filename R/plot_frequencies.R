#' plot genome-wide amplification and deletion frequencies
#'
#' @param cnr a cnr bundle
#'
#' @param genes logical, if TRUE use AmpFQ from gene data instead of bin, Default FALSE
#'
#' @param cols colors for amplifications and deletions, default to red and blue
#'
#' @param xlab x-axis label, passed to plot
#'
#' @param ylab y-axis label, passed to plot
#'
#' @param ylim sets coordinate space for y-axis, defaults to -1 to 1.  The
#' text will show absolute values for the deletion frequenies
#'
#' @param ... other arguments to plot
#' 
#' @return
#' Returns a bin alteration frequency plot with amplifications on the top axis,
#' and deletions on the bottom axis.
#' 
#' @examples
#'
#' data(cnr)
#' cnr <- get_alteration_frequencies(cnr)
#'
#' plot_frequencies(cnr)
#' plot_frequencies(cnr, xaxs = "i", bty = "n")
#'
#' @importFrom graphics plot abline lines axis
#'
#' 
#' @export
plot_frequencies <- function(cnr,
                             genes =  FALSE,
                             cols = c("#D40000", "#318CE7"),
                             xlab = "Genome", ylab = "Alteration Frequency",
                             ylim = c(-1, 1),...) {

    if(!genes) {
        AmpFQ <- cnr$chromInfo$AmpFQ
        delFQ <- cnr$chromInfo$delFQ *-1                                
    }
    if(genes) {
        AmpFQ <- cnr$gene.index$AmpFQ
        delFQ <- cnr$gene.index$delFQ *-1
    }

    chrBreaks <- chr_breaks(cnr, genes =  genes)
    midChr <- mid_chr(cnr, genes =  genes)

    graphics::plot(AmpFQ, type = "n", col = cols[1],
                   xaxt = "n", yaxt = "n", , ylim = ylim,
                   xlab = xlab, ylab = ylab, ...)
    graphics::abline(h = seq(-1, 1, by = 0.5), col = "gray")
    graphics::lines(AmpFQ, type = "h", col = cols[1], ...)
    graphics::lines(delFQ, type = "h", col = cols[2], ...)
    graphics::abline(h = 0)

    
    abline(v = c(0, chrBreaks), lty = 3, lwd = 0.8)
    graphics::axis(1, at = midChr, labels = names(midChr))
    graphics::axis(2, at = seq(-1, 1, by = 0.5), labels = abs(seq(-1, 1, by = 0.5)))
} # end plot_frequencies

