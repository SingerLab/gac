rm(list = ls())

library(testthat)
library(gac)

## load example data
data(copynumbers, pheno, qc, chromInfo, gene.index)

cnr <- buildCNR(X = copynumbers, Y = pheno, qc =qc,
                chromInfo = chromInfo, gene.index = gene.index)

expect_true("X" %in% names(cnr))
expect_true("genes" %in% names(cnr))
expect_true("Y" %in% names(cnr))
expect_true("qc" %in% names(cnr))
expect_true("chromInfo" %in% names(cnr))
expect_true("gene.index" %in% names(cnr))
expect_true("cells" %in% names(cnr))
expect_true("bulk" %in% names(cnr))

## check number of rows and columns throughout cnr object
expect_equal(length(cnr), 8)

expect_true(all(names(cnr) %in% c("X", "genes", "Y", "qc",
                                  "chromInfo", "gene.index", "cells", "bulk")))

## visualize genome-wide
h1 <- HeatmapCNR(cnr)
expect_true(all.equal(dim(h1@matrix), dim(cnr$X)))

## visualize genes of interest
h2 <- HeatmapCNR(cnr, what = "genes", which.genes = c("CDK4", "MDM2"))

## ADD cells
n.cells <- nrow(cnr$Y)
new.cells <- paste0("cell", n.cells+1:2)

newX <- data.frame(cbind(rep(c(5,2), c(3000, 2000)),
                         rep(c(2,4),  c(3000, 2000))))
names(newX) <- new.cells
head(newX)
tail(newX)

## creating new phenotypes
newY <- head(cnr$Y, n = 2)
newY$cellID <- new.cells
rownames(newY) <- newY$cellID

newY[, c(6:9)] <- newY[,c(6,9)]+3
head(newY)
tail(newY)

## creating new QC
newQC <- head(cnr$qc, 2)
newQC$cellID <- paste0("cell", 13:14)
rownames(newQC) <- newQC$cellID

newQC[,2:4] <- newQC[,2:4] + 2
head(newQC)
tail(newQC)

## add cells
cnr <- addCells(cnr, newX = newX, newY = newY, newqc = newQC)

## sapply(cnr, dim)

## expect_equal(length(cnr), 9)
expect_equal(ncol(cnr$X), length(cnr$cells))
expect_equal(nrow(cnr$Y), length(cnr$cells))
expect_equal(nrow(cnr$qc), length(cnr$cells))
expect_equal(length(cnr), 9)

h3 <- HeatmapCNR(cnr)
expect_equal(ncol(h3@matrix), n.cells +2)

## remove cells
( excl.cells <- rownames(cnr$qc)[cnr$qc$qc.status == "FAIL"] )
cnr <- excludeCells(cnr, excl = excl.cells)
sapply(cnr, dim)
n.cells2 <- ncol(cnr$X)

expect_equal(any(excl.cells %in% names(cnr$X)), FALSE)
expect_equal(any(excl.cells %in% rownames(cnr$Y)), FALSE)
expect_equal(any(excl.cells %in% rownames(cnr$qc)), FALSE)
expect_equal(any(excl.cells %in% cnr$cells), FALSE)

## expect_true(all(cnr$qc$qc.status == "PASS"))

h4 <- HeatmapCNR(cnr)

expect_equal(ncol(h4@matrix), n.cells2)

## keep cells
( keep.cells <- colnames(cnr$X)[c(1:8)] )
cnr <- keepCells(cnr, keep = keep.cells)
sapply(cnr, dim)

## addPheno
rand3 <- data.frame(cellID = cnr$Y$cellID,
                    rand3 = rnorm(nrow(cnr$Y), mean = 2, sd = 1))

cnr <- addPheno(cnr, df = rand3, by = "cellID", sort = FALSE)
expect_true(ncol(cnr$Y) == 10)

expect_true("rand3" %in% colnames(cnr$Y))

## add QC
mapd <- data.frame(t(apply(cnr$X, 2, mapd)))
mapd <- data.frame(cellID = rownames(mapd), mapd)

cnr <- addQC(cnr, df = mapd, by = "cellID", sort = FALSE)

expect_equal(ncol(cnr$qc), 9)
expect_true("mapd" %in% names(cnr$qc))
expect_true("mapd.sd" %in% names(cnr$qc))
expect_true("mapd.cv" %in% names(cnr$qc))


## addInfo

fakePval <- data.frame(pval = runif(5000))
cnr <- addInfo(cnr, df = fakePval)

expect_true("pval" %in% names(cnr$chromInfo))

