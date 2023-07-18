#' plot cluster optimization
#'
#' @param optC output of optClust()
#'
#' @param type type of plot e.g point, line, both (p, l, b)
#' 
#' @param ... additional arguments passed to plot and points
#'
#' @examples
#'
#' data(cnr)
#'
#' cnr <- phylo_cnr(cnr)
#'
#' ( mopc <- optClust(cnr, opt.range = seq(0.01, 0.2, by = 0.005)) )
#'
#' plotCL(optC = mopc)
#' 
#' @importFrom graphics plot points legend
#' 
#' @export
plotCL <- function(optC, type = "b", ...) {
    graphics::plot(rownames(optC), optC[,1], pch = 1, ylim = c(0, max(optC)),
         type = type, ...)
    graphics::points(rownames(optC), optC[,2], pch = 2, type = type,  ...)
    graphics::legend("topright", legend = c("One-cell", "Multi-cell"), pch = 1:2)
}


#' plot kCC vs sK
#' plot kParameter, and stable K from kStats
#' 
#' @param cnr a cnr bundle
#' 
#' @param kCC kCC, when NULL, kCC is taken from `optK`
#' 
#' @param sK sK, when NULL kS is taken from optK
#' 
#' @param pch pch for all points, default is 21
#' 
#' @param bg background color for all points, default is "gray"
#' 
#' @param xaxs x axis parameter see \code{\link{par}}, default is "i"
#' 
#' @param yaxs y axis parameter see \code{\link{par}}, default is "i"
#'
#' @param bty boxtype \code{\link{par}}, default is "l"
#'
#' @param xpd allow margin text margin text \code{\link{par}}, default is TRUE
#'
#' @param xlab x axis label, default is "k Parameter (kCC)"
#' 
#' @param ylab y axis label, default is "k Stable (kS)"
#' 
#' @param segment.col segment line color, default is Scotish blue "#005EB8"
#' 
#' @param segment.spacing segment spacing before highighted point, default is -0.8
#' 
#' @param segment.lwd segment line width, default 1.8
#' 
#' @param highlight.col color of a kCC /kS point, default is black "#282828"
#' 
#' @param highlight.cex character expansion of a kCC /kS point, default is 1.1
#' 
#' @param highlight.pch plot character of a kCC /kS point, default is 16
#'
#' @param ... other parameters passed to plot
#'
#' @importFrom assertthat assert_that
#' 
#' @importFrom graphics segments
#'
#' @export
plot_sK <- function(cnr, kCC = NULL, sK = NULL,
                    pch = 21, bg = "gray",
                    xaxs = "i", yaxs = "i", 
                    bty = "l", xpd = TRUE,
                    xlab = "k Parameter (kCC)", ylab = "k Stable (kS)",
                    segment.col = "#005EB8", segment.spacing = -0.8,
                    segment.lwd = 1.8,
                    highlight.col = "#282828",
                    highlight.cex = 1.1, highlight.pch = 16, ...) {

    assertthat::assert_that(!is.null(cnr$ccp), msg = "ccp not found, please run run_consensus_cluster")
    assertthat::assert_that(!is.null(cnr$kStats), msg = "kStats not found, please run kSpectral")
    assertthat::assert_that(!is.null(cnr$optK), msg = "optK not found, please run setKcc")
    if(is.null(kCC)) {
        kCC <- cnr$optK["kCC"]
    }
    if(is.null(sK)) {
        sK <- cnr$optK["sK"]
    }

    max.kcc <- max(cnr$kStats$kCC)
    max.sk <- max(cnr$kStats$k)
    
    plot(k ~ kCC, data = cnr$kStats, type = "l",
         xlim = c(0, max.kcc), ylim = c(0, max.sk),
         xaxs = xaxs, yaxs = yaxs, xpd = xpd,
         xlab = xlab, ylab = ylab,
         bty = bty,
         ... )
    
    graphics::segments(x0 = c(kCC, 0),
             x1 = c(kCC, kCC + segment.spacing),
             y0 = c(0, sK),
             y1 = c(sK + segment.spacing, sK),
             lwd = segment.lwd,
             col = segment.col,
             xpd = xpd)

    points(cnr$kStats[,c("kCC", "k")], 
           pch = pch, bg = bg, xpd = xpd)

    points(x = kCC, y = sK,
           col = highlight.col,
           cex = highlight.cex,
           pch = highlight.pch,
           xpd = xpd)

} #end plot_sK

