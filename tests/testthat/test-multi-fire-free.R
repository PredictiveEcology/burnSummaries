## Multi mode when no replicate had a fire. A fire-free replicate writes no parquet
## (fireregimetools::write_burn_parquet() skips an empty table), so multi mode has no schema to read.
## Through 1.1.0.9000 it built a table with no columns: burnSummaries_fireSizes_allReps.csv was
## written empty ("Input has no columns") and plotFun() stopped with "Object 'simArea' not found".

test_that("multi mode with every replicate fire-free keeps the fire-size columns", {
  mods <- writeLandMineStandIn(withr::local_tempdir())

  ## the columns and types multi mode reads back from a replicate with fires
  withFire <- file.path(testPaths$outputPath, "multi-withFire")
  runSingleRep(file.path(withFire, "rep01"), rep = 1L, fires = TRUE, mods = mods)
  ref <- data.table::as.data.table(dplyr::collect(
    fireregimetools::open_burn_dataset(file.path(withFire, "rep01", "burnSummaries_fireSizes"))
  ))
  expect_gt(nrow(ref), 0L)

  allFree <- file.path(testPaths$outputPath, "multi-fireFree")
  for (i in 1:2) {
    runSingleRep(file.path(allFree, sprintf("rep%02d", i)), rep = i, fires = FALSE, mods = mods)
  }
  expect_null(fireregimetools::open_burn_dataset(
    file.path(allFree, sprintf("rep%02d", 1:2), "burnSummaries_fireSizes")
  ))

  expect_no_warning(sim <- runMulti(allFree, reps = 1:2), message = "Input has no columns")

  expect_identical(nrow(sim$fireSizes), 0L)
  expect_identical(lapply(sim$fireSizes, class), lapply(ref, class))
  expect_identical(
    readLines(file.path(allFree, "burnSummaries_fireSizes_allReps.csv")),
    paste(names(ref), collapse = ",")
  )

  ## the burn maps are still drawn; the fire-size histogram, with no fires, is not
  figs <- basename(list.files(allFree, "\\.png$", recursive = TRUE))
  expect_true("cumulative_burn_maps.png" %in% figs)
  expect_false("simulated_number_size_fires.png" %in% figs)
})

test_that("multi mode writes the fire regime and fire per block, observed years limited by historicYears", {
  mods <- writeLandMineStandIn(withr::local_tempdir())
  out <- file.path(testPaths$outputPath, "multi-regime")
  runSingleRep(file.path(out, "rep01"), rep = 1L, fires = TRUE, mods = mods)
  runSingleRep(file.path(out, "rep02"), rep = 2L, fires = FALSE, mods = mods)
  sim <- runMulti(out, reps = 1:2, historicYears = c(1991L, 2000L))
  fr <- data.table::fread(file.path(out, "burnSummaries_fireRegime.csv"))
  expect_identical(fr$source, c("simulated", "simulated", "observed"))
  expect_identical(fr$nFires[1:2], c(4L, 0L)) ## LandMine stand-in: one fire each year 1-4
  ## the one observed polygon covers the 3 flammable cells once in the 10-year record
  pixHa <- prod(terra::res(testGrid())) / 1e4
  expect_equal(fr$areaBurnedHaPerYear[3], 3 * pixHa / 10)
  expect_identical(c(fr$firstYear[3], fr$lastYear[3]), c(1991L, 2000L))
  fb <- data.table::fread(file.path(out, "burnSummaries_fireBlocks.csv"))
  expect_identical(nrow(fb), 2L * 4L) ## 2 reps x 4 one-year blocks
  expect_true("fire_per_block.png" %in% basename(list.files(out, "\\.png$", recursive = TRUE)))
})

test_that("regimePeriod limits the simulated regime rows to those years", {
  mods <- writeLandMineStandIn(withr::local_tempdir())
  out <- file.path(testPaths$outputPath, "multi-regimePeriod")
  runSingleRep(file.path(out, "rep01"), rep = 1L, fires = TRUE, mods = mods)
  sim <- runMulti(out, reps = 1L, historicYears = c(1991L, 2000L), regimePeriod = c(2L, 4L))
  fr <- data.table::fread(file.path(out, "burnSummaries_fireRegime.csv"))
  expect_identical(fr$nFires[1], 2L) ## one fire each year 1-4; years 3 and 4 are in (2, 4]
  expect_identical(c(fr$firstYear[1], fr$lastYear[1]), c(2L, 4L))
})
