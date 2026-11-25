## LandMine (>= 1.0.13) records a fire that stalls in the second phase of its reburns as several rows of
## `fireSizes` -- what burned, then the shortfall ignited again from new start cells -- sharing the fire's
## `fireID` and original `targetSize`. `fireID` is issued per burn year. Through 1.1.0.9009 the fire
## regime and fire per block tables counted rows, so such a fire counted as several fires and the
## largest fire was its largest piece.

## rep 1, year 5: fire 1 in three pieces (6 + 3 + 1 of a target of 10) and fire 2 in one;
## rep 1, year 6: another fire 1; rep 2: two scfm-style rows with no fireID, one fire each.
piecesTable <- function() {
  data.frame(
    simArea = "test",
    rep = c(1L, 1L, 1L, 1L, 1L, 2L, 2L),
    year = c("5", "5", "5", "5", "6", "5", "6"),
    fireID = c(1L, 1L, 1L, 2L, 1L, NA, NA),
    attempt = c(3L, 4L, 4L, 1L, 1L, NA, NA),
    targetSize = c(10, 10, 10, 4, 7, NA, NA),
    simSize = c(6L, 3L, 1L, 4L, 7L, 2L, 5L),
    expSize = c(6L, 3L, 1L, 4L, 7L, NA, NA)
  )
}

test_that("a fire burned in pieces is one fire, its size the sum of its pieces", {
  f <- as.data.frame(firesFromPieces(piecesTable()))
  expect_identical(nrow(f), 5L)
  fire1 <- f[f$rep == 1L & f$year == "5" & f$fireID %in% 1L, ]
  expect_equal(fire1$simSize, 10)
  expect_equal(fire1$targetSize, 10)
  expect_identical(fire1$nPieces, 3L)
  ## the same fireID in another year is another fire
  expect_equal(f$simSize[f$rep == 1L & f$year == "6"], 7)
  ## a row with no fireID is one fire
  expect_equal(sort(f$simSize[f$rep == 2L]), c(2, 5))
  expect_true(all(f$nPieces[f$rep == 2L] == 1L))
})

test_that("the fire regime and fire per block count fires, not pieces", {
  f <- firesFromPieces(piecesTable())
  r <- simulatedFireRegime(f, reps = 1:2, simTimes = c(0, 10), pixelHa = 1)
  expect_identical(r$nFires, c(3L, 2L))
  expect_equal(r$largestFireHa, c(10, 5)) ## the whole fire, not its largest piece (6)
  ## area burned is the same counted either way
  expect_equal(r$areaBurnedHaPerYear, c(21, 7) / 10)
  b <- fireBlockTable(f, reps = 1:2, simTimes = c(0, 10), blockYears = 10, pixelHa = 1)
  expect_equal(b$nFires, c(3, 2))
  expect_equal(b$areaBurnedHa, c(21, 7))
})

test_that("a fireSizes table without the fire-identity columns has one fire per row", {
  fs <- piecesTable()[, c("simArea", "rep", "year", "simSize", "expSize")]
  f <- firesFromPieces(fs)
  expect_identical(nrow(f), nrow(fs))
  expect_true(all(is.na(f$targetSize)))
})

test_that("grouping the pieces leaves the table it was given unchanged", {
  fs <- data.table::as.data.table(piecesTable())
  before <- data.table::copy(fs)
  invisible(firesFromPieces(fs))
  expect_identical(fs, before)
})

test_that("an empty fireSizes table gives an empty table of fires", {
  fs <- data.table::data.table(
    simArea = character(0), rep = integer(0), year = character(0), simSize = integer(0),
    expSize = integer(0), fireID = integer(0), attempt = integer(0), targetSize = numeric(0)
  )
  f <- firesFromPieces(fs)
  expect_identical(nrow(f), 0L)
  r <- simulatedFireRegime(f, reps = 1L, simTimes = c(0, 10), pixelHa = 1)
  expect_identical(r$nFires, 0L)
})
