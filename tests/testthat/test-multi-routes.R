## Multi mode reads each replicate's fire sizes either from the CSV files `outputsDF` lists (as LandWeb runs it)
## or from the parquet partition in each replicate directory. Both routes, and the order the replicates are
## given in, must give the same tables; and a replicate's fire sizes count for the replicate whose directory
## they are in.

## The fire regime and fire per block tables multi mode wrote under `out`, in a fixed row order.
multiTables <- function(out) {
  lapply(c(regime = "fireRegime", blocks = "fireBlocks"), function(x) {
    d <- as.data.frame(data.table::fread(file.path(out, paste0("burnSummaries_", x, ".csv"))))
    d <- d[do.call(order, unname(d[intersect(c("source", "rep", "block"), names(d))])), , drop = FALSE]
    rownames(d) <- NULL
    d
  })
}

test_that("the CSV and parquet routes, and the order of the replicates, give the same tables", {
  mods <- writeLandMineStandIn(withr::local_tempdir())
  out <- file.path(testPaths$outputPath, "multi-routes")
  runSingleRep(file.path(out, "rep01"), rep = 1L, fires = TRUE, mods = mods)
  runSingleRep(file.path(out, "rep02"), rep = 2L, fires = FALSE, mods = mods)

  ## the stand-in's fire sizes and burn maps describe the same fires
  expect_no_warning(runMulti(out, reps = 1:2, plots = NA_character_), class = "burnSummariesBurnTotals")
  parquet <- multiTables(out)
  runMulti(out, reps = 1:2, outputsDF = repOutputs(out, 1:2), plots = NA_character_)
  csv <- multiTables(out)
  runMulti(out, reps = c(2L, 1L), plots = NA_character_)
  reversed <- multiTables(out)

  expect_identical(csv, parquet)
  expect_identical(reversed, parquet)
  sim <- parquet$regime[parquet$regime$source == "simulated", ]
  expect_identical(sim$nFires, c(4L, 0L)) ## one fire a year in replicate 1, none in 2
})

test_that("fire sizes tagged with another replicate's number stop multi mode", {
  mods <- writeLandMineStandIn(withr::local_tempdir())
  out <- file.path(testPaths$outputPath, "multi-wrongRep")
  ## a single-mode run in rep01 that was told it was replicate 2
  runSingleRep(file.path(out, "rep01"), rep = 2L, fires = TRUE, mods = mods)
  expect_error(runMulti(out, reps = 1L, plots = NA_character_), "rows for replicate 2 in the outputs of replicate 1")
})

test_that("a burn map that does not match the fire sizes warns, naming the replicate", {
  mods <- writeLandMineStandIn(withr::local_tempdir())
  out <- file.path(testPaths$outputPath, "multi-burnTotals")
  runSingleRep(file.path(out, "rep01"), rep = 1L, fires = TRUE, mods = mods)
  f <- file.path(out, "rep01", "burnMap_year4.tif")
  b <- terra::rast(f)
  v <- terra::values(b, mat = FALSE)
  v[2] <- 3 ## burns the fire sizes do not have
  terra::writeRaster(terra::rast(b, vals = v), f, overwrite = TRUE)
  expect_warning(
    runMulti(out, reps = 1L, plots = NA_character_),
    "1 \\(4 from fire sizes, 7 on the burn map\\)",
    class = "burnSummariesBurnTotals"
  )
})
