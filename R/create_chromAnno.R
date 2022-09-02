#' create chromosome annotations for custom heatmaps
#'
#' @param cnr a cnr bundle
#'
#' @param side side of the annotation. Options are "left" default,
#' "right", "top", and "bottom"
#'
#' @param ... additional prameters passed to HeatmapAnnotation
#'
#' @return
#' A HeatmapAnnotation containing a chromosome map indicate chromsome
#' boundaries on heatmaps
#'
#' @examples
#' data(cnr)
#'
#' chrAnnoLeft <- create_chromosome_annotation(cnr)
#' chrAnnoTop <- create_chromosome_annotation(cnr, side = "top")
#'
#' @export
create_chromosome_annotation <- function(cnr, side = "left", ...) {
    if(side == "left") {
        chA <- create_chromosome_annotation_left(cnr, ...)
    }

    if(side == "right") {
        chA <- create_chromosome_annotation_right(cnr, ...)
    }

    if(side == "top") {
        chA <- create_chromosome_annotation_top(cnr, ...)
    }

    if(side == "bottom") {
        chA <- create_chromosome_annotation_bottom(cnr, ...)
    }

    return(chA)

}


#' create chromosome annotations for custom heatmaps
#'
#' @param cnr a cnr bundle
#'
#' @param ... additional prameters passed to HeatmapAnnotation
#' 
#' @importFrom ComplexHeatmap rowAnnotation anno_mark
#' @importFrom grid gpar
create_chromosome_annotation_left <- function(cnr, ...) {
    cf <- factor(cnr$chromInfo$bin.chrom)
    grs <- c("#404040", "#BABABA")
    rp <- ceiling(length(unique(cf))/2)
    chl <- rep(grs, rp)
    chl <- chl[1:length(unique(cf))]
    names(chl) <- unique(cf)

    chrBreaks <- cumsum(table(cnr$chromInfo$bin.chrom))
    if (length(chrBreaks) == 1) {
        midChr <- floor(chrBreaks/2)
    }
    else {
        midChr <-
            chrBreaks - floor((chrBreaks - c(1, chrBreaks[1:(length(chrBreaks) - 
                                                             1)]))/2)
    }
    chrAnno <- rowAnnotation(
        labs = anno_mark(at = midChr, 
                         labels = unique(cf),
                         side = "left",
                         labels_gp = grid::gpar(fontsize = 10)), 
        chr = cf, col = list(chr = chl),
        show_annotation_name = FALSE,
        show_legend = FALSE)


    return(chrAnno)
}

#' create chromosome annotations for custom heatmaps
#'
#' @param cnr a cnr bundle
#' @param labels_gp graphic parameters from \link[grid]{gpar}, default fontsize = 10
#' @param labels_rot label rotation, default 90
#' @param ... additional prameters passed to HeatmapAnnotation
#' 
#' @importFrom ComplexHeatmap rowAnnotation anno_mark
#' @importFrom grid gpar
create_chromosome_annotation_top <- function(cnr,
                                             labels_gp = grid::gpar(fontsize = 10),
                                             labels_rot = 90, ...) {
    cf <- factor(cnr$chromInfo$bin.chrom)
    grs <- c("#404040", "#BABABA")
    rp <- ceiling(length(unique(cf))/2)
    chl <- rep(grs, rp)
    chl <- chl[1:length(unique(cf))]
    names(chl) <- unique(cf)

    chrBreaks <- cumsum(table(cnr$chromInfo$bin.chrom))
    if (length(chrBreaks) == 1) {
        midChr <- floor(chrBreaks/2)
    }
    else {
        midChr <-
            chrBreaks - floor((chrBreaks - c(1, chrBreaks[1:(length(chrBreaks) - 
                                                             1)]))/2)
    }
    chrAnno <- HeatmapAnnotation(
        labs = anno_mark(at = midChr, 
                         labels = unique(cf),
                         side = "top",
                         labels_gp = labels_gp,
                         labels_rot = labels_rot, ...), 
        chr = cf, col = list(chr = chl),
        show_annotation_name = FALSE,
        show_legend = FALSE)

    return(chrAnno)
}

#' create chromosome annotations for custom heatmaps
#'
#' @param cnr a cnr bundle
#' @param labels_gp graphic parameters from \link[grid]{gpar}, default fontsize = 10
#' @param labels_rot label rotation, default 90
#' @param ... additional prameters passed to HeatmapAnnotation
#' 
#' @importFrom ComplexHeatmap rowAnnotation anno_mark
#' @importFrom grid gpar
create_chromosome_annotation_bottom <- function(cnr,
                                             labels_gp = grid::gpar(fontsize = 10),
                                             labels_rot = 90, ...) {
    cf <- factor(cnr$chromInfo$bin.chrom)
    grs <- c("#404040", "#BABABA")
    rp <- ceiling(length(unique(cf))/2)
    chl <- rep(grs, rp)
    chl <- chl[1:length(unique(cf))]
    names(chl) <- unique(cf)

    chrBreaks <- cumsum(table(cnr$chromInfo$bin.chrom))
    if (length(chrBreaks) == 1) {
        midChr <- floor(chrBreaks/2)
    }
    else {
        midChr <-
            chrBreaks - floor((chrBreaks - c(1, chrBreaks[1:(length(chrBreaks) - 
                                                             1)]))/2)
    }
    chrAnno <- HeatmapAnnotation(
        chr = cf, col = list(chr = chl),
        labs = anno_mark(at = midChr, 
                         labels = unique(cf),
                         side = "bottom",
                         labels_gp = labels_gp,
                         labels_rot = labels_rot, ...), 
        show_annotation_name = FALSE,
        show_legend = FALSE)

    return(chrAnno)
}

#' create chromosome annotations for custom heatmaps
#'
#' @param cnr a cnr bundle
#' @param ... additional prameters passed to HeatmapAnnotation
#' 
#' @importFrom ComplexHeatmap rowAnnotation anno_mark
#' @importFrom grid gpar
create_chromosome_annotation_right <- function(cnr, ...) {
    cf <- factor(cnr$chromInfo$bin.chrom)
    grs <- c("#404040", "#BABABA")
    rp <- ceiling(length(unique(cf))/2)
    chl <- rep(grs, rp)
    chl <- chl[1:length(unique(cf))]
    names(chl) <- unique(cf)

    chrBreaks <- cumsum(table(cnr$chromInfo$bin.chrom))
    if (length(chrBreaks) == 1) {
        midChr <- floor(chrBreaks/2)
    }
    else {
        midChr <-
            chrBreaks - floor((chrBreaks - c(1, chrBreaks[1:(length(chrBreaks) - 
                                                             1)]))/2)
    }

    chrAnno <- rowAnnotation(
        chr = cf, col = list(chr = chl),
        labs = anno_mark(at = midChr, 
                         labels = unique(cf),
                         side = "right",
                         labels_gp = grid::gpar(fontsize = 10),
                         ...), 
        show_legend = FALSE)
    
    return(chrAnno)
}
