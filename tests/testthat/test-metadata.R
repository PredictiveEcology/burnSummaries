## The module's metadata is its public contract: a project using this module binds to these
## object names and classes, and `reqdPkgs` states what it needs to run at all. These are
## CHARACTERIZATION tests -- they pin today's contract so a change to it has to be deliberate,
## rather than describing behaviour that did not exist before.
##
## GENERATED from the module's live metadata, then reviewed. When a change is intended, update
## this file in the same commit and bump the module version to match: removed, renamed or
## retyped is a MAJOR bump.

test_that("module metadata parses", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_type(md, "list")
  expect_identical(md$name, moduleName)
})

test_that("inputs are the expected names and classes", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  inputs <- stats::setNames(md$inputObjects$objectClass, md$inputObjects$objectName)
  inputs <- inputs[!is.na(names(inputs))]
  expect_identical(
    inputs[order(tolower(names(inputs)))],
    c(
      "burnMap"                        = "SpatRaster",
      "burnSummary"                    = "data.table",
      "fireSizes"                      = "list",
      "flammableMap"                   = "SpatRaster",
      "nonForest_timeSinceDisturbance" = "SpatRaster",
      "rstCurrentBurn"                 = "SpatRaster",
      "rstTimeSinceFire"               = "SpatRaster"
)
  )
})

test_that("outputs are the expected names and classes", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  outputs <- stats::setNames(md$outputObjects$objectClass, md$outputObjects$objectName)
  outputs <- outputs[!is.na(names(outputs))]
  expect_identical(
    outputs[order(tolower(names(outputs)))],
    c(
      "fireSizes"        = "data.table",
      "rstTimeSinceFire" = "SpatRaster"
)
  )
})

test_that("parameters are the expected names", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_identical(
    sort(md$parameters$paramName),
    c(
      ".plotInitialTime", ".plotInterval", ".plots", ".saveInitialTime",
      ".saveInterval", ".seed", ".studyAreaName", ".useCache", "dataYear",
      "fireTimestep", "mode", "reps", "simOutputPath", "simTimes", "summaryInterval",
      "summaryPeriod"
)
  )
})

