test_that("addPheno appends a new column to Y", {
  data(cnr)
  n_cols <- ncol(cnr$Y)
  set.seed(42)
  new_col <- data.frame(cellID = cnr$Y$cellID, rand3 = rnorm(nrow(cnr$Y)))

  cnr2 <- addPheno(cnr, df = new_col, by = "cellID", sort = FALSE)

  expect_equal(ncol(cnr2$Y), n_cols + 1)
  expect_true("rand3" %in% colnames(cnr2$Y))
})

test_that("addQC appends mapd columns to qc", {
  data(cnr)
  mapd_df <- data.frame(t(apply(cnr$X, 2, mapd)))
  mapd_df <- data.frame(cellID = rownames(mapd_df), mapd_df)

  cnr2 <- addQC(cnr, df = mapd_df, by = "cellID", sort = FALSE)

  expect_true("mapd" %in% names(cnr2$qc))
  expect_true("mapd.sd" %in% names(cnr2$qc))
  expect_true("mapd.cv" %in% names(cnr2$qc))
})

test_that("addInfo appends columns to chromInfo", {
  data(cnr)
  fake_pval <- data.frame(pval = runif(nrow(cnr$chromInfo)))

  cnr2 <- addInfo(cnr, df = fake_pval)

  expect_true("pval" %in% names(cnr2$chromInfo))
})
