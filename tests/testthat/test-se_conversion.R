test_that("conversion to DGEList works", {
  library(SummarizedExperiment)
  se = SummarizedExperiment(assays = list('counts' = matrix(0, 3, 3)))

  expect_equal(dim(asDGEList(se)), dim(se))
  expect_error(dim(asDGEList(se, 'qwerty')))
  #expect error when no assay data present
  expect_error(dim(asDGEList(SummarizedExperiment())))
  expect_error(dim(asDGEList(SummarizedExperiment(), 1)))
})
