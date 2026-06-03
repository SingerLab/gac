test_that("buildCNR produces all required slots", {
  data(copynumbers, pheno, qc, chromInfo, grch37.genes.5k)
  cnr <- buildCNR(X = copynumbers, Y = pheno, qc = qc,
                  chromInfo = chromInfo, gene.index = grch37.genes.5k)

  expect_equal(length(cnr), 8)
  expect_true(all(c("X", "genes", "Y", "qc", "chromInfo", "gene.index", "cells", "bulk") %in% names(cnr)))
})

test_that("buildCNR keeps all matrices synchronized", {
  data(copynumbers, pheno, qc, chromInfo, grch37.genes.5k)
  cnr <- buildCNR(X = copynumbers, Y = pheno, qc = qc,
                  chromInfo = chromInfo, gene.index = grch37.genes.5k)

  expect_equal(ncol(cnr$X), nrow(cnr$Y))
  expect_equal(ncol(cnr$X), nrow(cnr$qc))
  expect_equal(ncol(cnr$X), nrow(cnr$genes))
  expect_equal(ncol(cnr$X), length(cnr$cells))
})

test_that("HeatmapCNR on bins matches X dimensions", {
  data(cnr)
  h <- HeatmapCNR(cnr)
  expect_equal(dim(h@matrix), dim(cnr$X))
})

test_that("HeatmapCNR on genes runs without error", {
  data(cnr)
  expect_no_error(HeatmapCNR(cnr, what = "genes", which.genes = c("CDK4", "MDM2")))
})
