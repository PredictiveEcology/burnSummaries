test_that("repFiles() matches each replicate directory and file name exactly", {
  files <- file.path("out", c("rep1", "rep10", "rep2"),
                     rep(c("burnMap_year1000.tif", "burnSummaries_fireSizes.csv"), each = 3))
  ## other years' maps, and files in a replicate's subfolders (as outputs(sim) lists them)
  files <- c(files, "out/rep1/burnMap_year0500.tif", "out/rep1/figures/burnMap_year1000.tif")
  repDirs <- dirnamesFromSet(files, 1:10, leafOnly = FALSE)
  expect_identical(repDirs, file.path("out", c("rep1", "rep2", "rep10")))

  ## rep1 does not pick up rep10's files, other years' maps or files in subfolders
  expect_identical(repFiles(repDirs, "burnMap_year1000.tif", files),
                   file.path(repDirs, "burnMap_year1000.tif"))
  expect_identical(repFiles(repDirs, "burnSummaries_fireSizes.csv", files),
                   file.path(repDirs, "burnSummaries_fireSizes.csv"))
})

test_that("repFiles() stops when a replicate's file is not in outputsDF", {
  files <- file.path("out", c("rep1", "rep2"), "burnMap_year1000.tif")[1]
  expect_error(repFiles(file.path("out", c("rep1", "rep2")), "burnMap_year1000.tif", files),
               "expected one burnMap_year1000.tif in out/rep2")
})

test_that("repFiles() without outputsDF builds the paths", {
  expect_identical(repFiles(file.path("out", c("rep01", "rep02")), "burnMap_year1000.tif"),
                   file.path("out", c("rep01", "rep02"), "burnMap_year1000.tif"))
})
