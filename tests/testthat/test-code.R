test_that("the module code contains no browser() call", {
  parsed <- parse(file.path(moduleRoot, paste0(moduleName, ".R")), keep.source = TRUE)
  expect_false("browser" %in% unlist(lapply(parsed, all.names)))
})
