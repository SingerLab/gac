#' Create chromosome annotations for custom heatmaps
#'
#' @param cnr a cnr bundle
#'
#' @param side side of the annotation. Options are "left" default,
#' "right", "top", and "bottom"
#'
#' @param ... additional prameters passed to \code{\link[ComplexHeatmap]{HeatmapAnnotation}}
#'
#' @return
#' A HeatmapAnnotation containing a chromosome map indicates chromosome
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


#' Create chromosome annotations for custom heatmaps
#'
#' @param cnr a cnr bundle
#'
#' @param labels_gp graphic parameters from \link[grid]{gpar}, default fontsize = 10
#' 
#' @param labels_rot label rotation, default 0
#'
#' @param ... additional parameters passed to HeatmapAnnotation
#' 
#' @importFrom ComplexHeatmap rowAnnotation anno_mark
#' @importFrom grid gpar
create_chromosome_annotation_left <- function(cnr,
                                              labels_gp = grid::gpar(fontsize = 10),
                                              labels_rot = 0, ...) {

    if(is.factor(cnr$chromInfo$bin.chrom)) {
        cf <- droplevels(cnr$chromInfo$bin.chrom)
    } else {
        cf <- factor(cnr$chromInfo$bin.chrom)
    }

    midChr <- mid_chr(cnr)
    chl <- chr_colors(cnr)
    
    chrAnno <- ComplexHeatmap::rowAnnotation(
        labs = ComplexHeatmap::anno_mark(at = midChr, 
                         labels = names(midChr),
                         side = "left",
                         labels_gp = labels_gp,
                         labels_rot = labels_rot), 
        chr = cf, col = list(chr = chl[names(midChr)]),
        show_annotation_name = FALSE,
        show_legend = FALSE,
        ...)

    return(chrAnno)
}

#' Create chromosome annotations for custom heatmaps
#'
#' @param cnr a cnr bundle
#' @param labels_gp graphic parameters from \link[grid]{gpar}, default fontsize = 10
#' @param labels_rot label rotation, default 0
#' @param ... additional parameters passed to HeatmapAnnotation
#' 
#' @importFrom ComplexHeatmap HeatmapAnnotation anno_mark
#' @importFrom grid gpar
create_chromosome_annotation_top <- function(cnr,
                                             labels_gp = grid::gpar(fontsize = 10),
                                             labels_rot = 0, ...) {
    if(is.factor(cnr$chromInfo$bin.chrom)) {
        cf <- droplevels(cnr$chromInfo$bin.chrom)
    } else {
        cf <- factor(cnr$chromInfo$bin.chrom)
    }
    
    midChr <- mid_chr(cnr)
    chl <- chr_colors(cnr)
    
    chrAnno <- ComplexHeatmap::HeatmapAnnotation(
        labs = ComplexHeatmap::anno_mark(at = midChr, 
                                         labels = names(midChr),
                         side = "top",
                         labels_gp = labels_gp,
                         labels_rot = labels_rot), 
        chr = cf, col = list(chr = chl[names(midChr)]),
        show_annotation_name = FALSE,
        show_legend = FALSE,
        ...)

    return(chrAnno)
}

#' Create chromosome annotations for custom heatmaps
#'
#' @param cnr a cnr bundle
#' @param labels_gp graphic parameters from \link[grid]{gpar}, default fontsize = 10
#' @param labels_rot label rotation, default 0
#' @param ... additional parameters passed to HeatmapAnnotation
#' 
#' @importFrom ComplexHeatmap HeatmapAnnotation anno_mark
#' @importFrom grid gpar
create_chromosome_annotation_bottom <- function(cnr,
                                             labels_gp = grid::gpar(fontsize = 10),
                                             labels_rot = 0, ...) {
    
    if(is.factor(cnr$chromInfo$bin.chrom)) {
        cf <- droplevels(cnr$chromInfo$bin.chrom)
    } else {
        cf <- factor(cnr$chromInfo$bin.chrom)
    }
    
    midChr <- mid_chr(cnr)
    chl <- chr_colors(cnr)
    
    chrAnno <- ComplexHeatmap::HeatmapAnnotation(
        chr = cf, col = list(chr = chl[names(midChr)]),
        labs = ComplexHeatmap::anno_mark(at = midChr, 
                         labels = names(midChr),
                         side = "bottom",
                         labels_gp = labels_gp,
                         labels_rot = labels_rot), 
        show_annotation_name = FALSE,
        show_legend = FALSE,
        ...)

    return(chrAnno)
}

#' Create chromosome annotations for custom heatmaps
#'
#' @param cnr a cnr bundle
#'
#' @param labels_gp graphic parameters from \link[grid]{gpar}, default fontsize = 10
#'
#' @param labels_rot label rotation, default 0
#'
#' @param ... additional parameters passed to HeatmapAnnotation
#' 
#' @importFrom ComplexHeatmap rowAnnotation anno_mark
#' @importFrom grid gpar
create_chromosome_annotation_right <- function(cnr,
                                               labels_gp = grid::gpar(fontsize = 10),
                                               labels_rot = 0, ...) {

    if(is.factor(cnr$chromInfo$bin.chrom)) {
        cf <- droplevels(cnr$chromInfo$bin.chrom)
    } else {
        cf <- factor(cnr$chromInfo$bin.chrom)
    }
    
    midChr <- mid_chr(cnr)
    chl <- chr_colors(cnr)
    
    chrAnno <- ComplexHeatmap::rowAnnotation(
        chr = cf, col = list(chr = chl[names(midChr)]),
        labs = ComplexHeatmap::anno_mark(at = midChr, 
                         labels = names(midChr),
                         side = "right",
                         labels_gp = labels_gp,
                         labels_rot = labels_rot),
        show_legend = FALSE,
        ...)
    
    return(chrAnno)
}


#' Estimate chromosome midpoint locations along a continuous genome
#'
#' @param cnr a cnr
#'
#' @param bin  weather to use bin or gene data, default is true
#'
#' @return
#' A named vector of chromosome midpoints. Useful for adding tick
#' marks in figures.  Midpoint is not the centromere location.
#' 
#' @export
mid_chr <- function(cnr, bin = TRUE) {
    brk <- chr_breaks(cnr, bin = bin)

    if (length(brk) == 1) {
        mid.pt <- floor(brk/2)
    } else {
        mid.pt <-
            brk - floor((brk - c(1, brk[1:(length(brk) - 1)]))/2)
    }
        return(mid.pt)
}


#' Estimate chromosome end locations along a continuous genome
#'
#' @param cnr a cnr
#'
#' @param bin  weather to use bin or gene data, default is true
#'
#' @return
#' A named vector of chromosome breaks locations in the data.
#' Useful when adding lines to separate chromosomes,  or
#' a background when highlighting a chromosome
#' 
#' @export
chr_breaks <- function(cnr, bin = TRUE) {

    if(is.factor(cnr$chromInfo$bin.chrom)) {
        cnr$chromInfo$bin.chrom <- droplevels(cnr$chromInfo$bin.chrom)
    } else {
        cnr$chromInfo$bin.chrom <- factor(cnr$chromInfo$bin.chrom)
    }

    if(bin) {
        brk <- cumsum(table(cnr$chromInfo$bin.chrom))
    } else {
        brk <- cumsum(table(cnr$gene.index$chrom))
    }

    return(brk)
}


#' Generate chromosome sidebar colors
#' @param cnr a cnr bundle
#'
#' @param col alternating chromosome colors, default is c("#404040", "#BABABA")
#' 
#' @param bin  weather to use bin or gene data, default is true
#'
#' @return
#' A named vector of default chromosome colors
#' @export
chr_colors <- function(cnr, col = c("#404040", "#BABABA"),
                       bin = TRUE) {

    if(is.factor(cnr$chromInfo$bin.chrom)) {
        cf <- droplevels(cnr$chromInfo$bin.chrom)
    } else {
        cf <- factor(cnr$chromInfo$bin.chrom)
    }
    
    rp <- ceiling(length(unique(cf))/2)

    chl <- rep(col, rp)
    chl <- chl[1:length(unique(cf))]
    names(chl) <- unique(cf)

    return(chl)
}

