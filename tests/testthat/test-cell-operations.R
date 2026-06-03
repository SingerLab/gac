test_that("addCells increases cell count across all slots", {
  data(cnr)
  n_start <- ncol(cnr$X)
  new.cells <- paste0("cell_new_", 1:2)

  newX <- data.frame(matrix(2L, nrow = nrow(cnr$X), ncol = 2))
  names(newX) <- new.cells

  newY <- head(cnr$Y, 2)
  newY$cellID <- new.cells
  rownames(newY) <- new.cells

  newQC <- head(cnr$qc, 2)
  newQC$cellID <- new.cells
  rownames(newQC) <- new.cells

  cnr2 <- addCells(cnr, newX = newX, newY = newY, newqc = newQC)

  expect_equal(ncol(cnr2$X), n_start + 2)
  expect_equal(nrow(cnr2$Y), n_start + 2)
  expect_equal(nrow(cnr2$qc), n_start + 2)
  expect_equal(length(cnr2$cells), n_start + 2)
})

test_that("excludeCells removes cells from all slots", {
  data(cnr)
  excl <- rownames(cnr$qc)[cnr$qc$qc.status == "FAIL"]
  n_start <- ncol(cnr$X)

  cnr2 <- excludeCells(cnr, excl = excl)

  expect_equal(ncol(cnr2$X), n_start - length(excl))
  expect_false(any(excl %in% colnames(cnr2$X)))
  expect_false(any(excl %in% rownames(cnr2$Y)))
  expect_false(any(excl %in% rownames(cnr2$qc)))
  expect_false(any(excl %in% cnr2$cells))
})

test_that("keepCells retains only specified cells", {
  data(cnr)
  keep <- colnames(cnr$X)[1:8]

  cnr2 <- keepCells(cnr, keep = keep)

  expect_equal(ncol(cnr2$X), 8)
  expect_equal(nrow(cnr2$Y), 8)
  expect_equal(nrow(cnr2$qc), 8)
  expect_equal(length(cnr2$cells), 8)
  expect_true(all(colnames(cnr2$X) %in% keep))
})
