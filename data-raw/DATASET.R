## code to prepare `DATASET` dataset goes here
## usethis::use_data(DATASET, overwrite = TRUE)

#' For convenience, GAC comes with several color palettes as part of the data.
#'  To use any of these run `data(segCol)` for example.  There are the different
#'  color palettes such as:
#' 
#'   * segCol for segmentation colors: loss and deletions in blue and yellow,
#'   respectively, 2 in white, gains and amplifications in red to gray scale
#' 
#'   * lowCol for log2 ratio colors: blue-to-red scale from -2 to +2
#' 
#'   * ampCol for amplification : yellow to red scale from 0 to Inf.
#'
#'   * fgaCol for Fraction of genome altered: from 0 to 1 light yellow to 1 in blue
#'
#'   * ploidyCol for ploidy : gradient of purples starting at 1.5
#'
#'   * ccClustCol : yellow to blue from 0 to 1 for co-clustering frequencies

devtools::load_all()

## fixed colors
segCol <- circlize::colorRamp2(breaks = c(0:5, 10, 20, 40, 80),
                     colors = c("#F2D200", "#318CE7", "#FFFFFF",
                                "#FFA89F", "#FF523F", "#D40000",
                                "#7F7F7F", "#636363", "#515151",
                                "#303030"))
lowCol <- circlize::colorRamp2(breaks = c(-2, 0, 2),
                     colors = c("#032ADD", "#FFFFFF", "#D40000"))
ampCol <- segCol

fgaCol <- circlize::colorRamp2(0:50/100, colourvalues::colour_values(0:50/100, "GnBu"))

ploidyCol <- circlize::colorRamp2(c(0, 2, 4),
                        colors = c("#F2F0F7", "#9E9AC8", "#54278F"))

## legends
legSeg <- Legend(at = c(0:5, 10, 20, 40, 80),
                 labels = c(0:5, 10, 20, 40, ">80"),
                 legend_gp = gpar(fill = colorRampPalette(
                                      c("#F2D200", "#318CE7", "#FFFFFF",
                                        "#FFA89F", "#FF523F", "#D40000",
                                        "#7F7F7F", "#636363", "#515151",
                                        "#303030"))(10)), title = "GeneCN")
legLog2R <- Legend(at = seq(-2, 2, by = 1), col_fun = lowCol, title = "Log2R")
legAmp <- legSeg
legFGA <- Legend(at = seq(0, .50, by = .10), col_fun = fgaCol, title = "FGA")
legPloidy <- Legend(at = c(0, 2, 4), col_fun = ploidyCol, title = "Ploidy")

usethis::use_data(segCol, lowCol, ampCol, fgaCol, ploidyCol,
                  legSeg, legLog2R, legAmp, legFGA, legPloidy,
                  overwrite = TRUE, compress = "xz")

## read data for cnr
X <- read.delim("data-raw/copynumbers.txt", as.is = TRUE)
Y <- read.delim("data-raw/pheno.txt", as.is = TRUE)
qc <- read.delim("data-raw/qc.txt", as.is = TRUE)
gx <- read.delim("data-raw/gene.index.txt", as.is = TRUE, na.strings = c(NA, ""))
gx <- gx[!is.na(gx$hgnc.symbol),]
gx <- gx[!(duplicated(gx$hgnc.symbol)),]
rownames(gx) <- gx$hgnc.symbol
ci <- read.delim("data-raw/chromInfo.txt", as.is = TRUE)

cnr <- buildCNR(X = X, Y = Y, qc = qc, exprs = NULL,
                gene.index = gx, chromInfo = ci)

HeatmapCNR(cnr, col = segCol)

dna <- buildCNR(X = X, Y = Y, qc = qc, exprs = NULL,
                gene.index = gx, chromInfo = ci, bulk = TRUE)

HeatmapCNR(dna)

HeatmapCNR(dna, col = lowCol)

usethis::use_data(cnr, dna, overwrite = TRUE, compress = "xz")

gene.index <- gx
chromInfo <- ci
copynumbers <- X

usethis::use_data(copynumbers, Y, qc, chromInfo, gene.index,
                  overwrite = TRUE, compress = "xz")


## currently in development -- use of package dm to speficy data model and
## the relationship between each table
## cnr2 <- dmCNR(X, Y, qc, ci, gx, bulk = FALSE)
