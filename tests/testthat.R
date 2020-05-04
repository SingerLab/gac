library(testthat)
library(toSignac)

## import example cnr
data(cnr)

## import colors
data(segCol)

## check number of rows and columns throughout cnr object
sapply(cnr, dim)

## visualize genome-wide
HeatmapCNR(cnr)


## visualize genes of interest
HeatmapCNR(cnr, what = "genes", which.genes = c("CDK4", "MDM2"))

## ADD cells

newX <- data.frame(cbind(rep(c(5,2), c(3000, 1999)),
                         rep(c(2,4),  c(3000, 1999))))
names(newX) <- paste0("cell", 13:14)
newX

## creating new phenotypes
newY <- head(cnr$Y, n = 2)
newY$cellID <- paste0("cell", 13:14)
rownames(newY) <- newY$cellID

newY[, c(6:9)] <- newY[,c(6,9)]+3
newY

## creating new QC
newQC <- head(cnr$qc, 2)
newQC$cellID <- paste0("cell", 13:14)
rownames(newQC) <- newQC$cellID

newQC[,2:4] <- newQC[,2:4] + 2
newQC

## add cells
cnr <- addCells(cnr, newX = newX, newY = newY, newqc = newQC)
sapply(cnr, dim)
HeatmapCNR(cnr)

## remove cells
( excl.cells <- rownames(cnr$qc)[cnr$qc$qc.status == "FAIL"] )
cnr <- excludeCells(cnr, excl = excl.cells)
sapply(cnr, dim)
HeatmapCNR(cnr)

## keep cells
( keep.cells <- colnames(cnr$X)[-c(11:12)] )
cnr <- keepCells(cnr, keep = keep.cells)
sapply(cnr, dim)

HeatmapCNR(cnr)


## addPheno
rand3 <- data.frame(cellID = cnr$Y$cellID, rand3 = rnorm(nrow(cnr$Y), mean = 2, sd = 1))

cnr <- addPheno(cnr, df = rand3)
sapply(cnr, dim)

## add QC
mapd <- data.frame(t(apply(cnr$X, 2, mapd)))
mapd <- data.frame(cellID = rownames(mapd), mapd)

cnr <- addQC(cnr, df = mapd)

## addInfo

fakePval <- data.frame(runif(4999))

cnr <- addInfo(cnr, df = fakePval)

