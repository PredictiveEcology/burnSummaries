# burnSummaries (development version)

## Reuse the cached stand-age input across replicates (`1.0.2.9007`)

* `.inputObjects` now downloads/reads the SCANFI stand-age source used to seed `rstTimeSinceFire`
  into `inputPath(sim)` (the shared inputs cache) instead of the per-replicate `outputPath(sim)`.
  Previously every replicate re-downloaded the ~5.2 GB SCANFI age file to its own output directory,
  so a multi-replicate mainSim launched many simultaneous large downloads; one dropped its
  connection and, with the headless no-retry guard, failed the whole run. Pointing at the shared
  inputs path reuses the already-cached file (no re-download).

## Robust downloads for the large NBAC/NFDB archives (`1.0.2.9006`)

* Downloading the ~1.2 GB NBAC composite exceeded R's default 60 s `download.file` timeout, silently
  truncating the zip (extraction then failed with `ZIP decompression failed`). Raise the timeout and
  download to a `.part` file that is renamed only on success, so a truncated/interrupted download is
  not mistaken for a complete one on a later run.

## Self-contained summary-output times (`1.0.2.9005`)

* Inline the summary-output-times calculation (`seq()` over the summary period) instead of calling
  `LandWebUtils::analysesOutputsTimes()`. burnSummaries is a generic module and does not declare
  LandWebUtils; the bare call only resolved when a LandWeb module (e.g. NRV_summary) was co-run and
  loaded it, so a standalone `mode = "multi"` run errored with `could not find function`.

## Observed fire perimeters from NBAC + NFDB backfill (`1.0.2.9004`)

* The historical (observed) cumulative burn map now uses **National Burned Area Composite (NBAC)**
  perimeters (satellite-derived, 1972-present) as the authoritative source, supplemented with
  **National Fire DataBase (NFDB)** polygons ONLY for years NBAC does not cover. Older NFDB
  perimeters are aerial sketches that overestimate burned area, so NBAC is preferred wherever it
  exists. Loading + harmonising (tolerant `YEAR`/`SIZE_HA` columns, clipped to the sim grid) is done
  via `fireregimetools::load_nbac_polys()` / `load_nfdb_polys()` (>= 0.1.0), replacing the former
  NFDB-only historical burn map.

## Fire-regime summaries via fireregimetools (`1.0.2.9003`)

* Adopt the shared, arrow-native `FOR-CAST/fireregimetools` package for the fire-size summaries.
  `create_fireSizes` now also publishes each replicate's fire-size table as a parquet partition
  (`fireregimetools::write_burn_parquet`); `multi`-mode `FireSummaries` reads all replicates as one
  lazy Arrow dataset (`fireregimetools::open_burn_dataset`) instead of `rbind`-ing the per-replicate
  CSVs into memory; and the fire-size distribution plots (`ggHistSim` / `ggHistExp`) are drawn by
  `fireregimetools::fire_size_histogram()` (count histogram with a median-log-size-per-bin overlay),
  replacing the bespoke `hist()` + `stat_summary_bin()` dual-axis code. The per-replicate and
  all-reps CSVs are still written. Adds `FOR-CAST/fireregimetools` to `reqdPkgs`.

*NEWS was not maintained between the initial `0.0.1` module and the current development
version (`1.0.2.9001`); this entry catches up the substantive changes over that window.*

## Module architecture and modes

* Extracted the burn post-processing code out of the main LandWeb repository into a
  self-contained module (2023-09); reached `1.0.0` alongside an scfm update and adoption of
  `terra` (2024-06).
* Added a two-phase single/multi mode. `single` mode saves the inputs that `multi` mode needs
  as its own outputs (avoiding unreliable whole-`simList` loads), incorporates the former
  `timeSinceFire` module, calls `registerOutputs()`, and switches serialization to `qs2`;
  `multi` mode loads from the saved files. Plotting moved from `rasterVis` to `ggplot2` +
  `tidyterra`, and spatial objects to `terra` vectors.
* Added custom per-object saving in `single` mode; expected/simulated expected-value
  calculation and plotting are skipped when target fire sizes are absent.

## Simulation-time alignment

* `analysesOutputsTimes` and the scheduled `save_single` events are now offset by `start(sim)`
  so summary times are simulation-start-relative.
* Added an `Init`-time guard that `summaryPeriod` falls within `start(sim)`/`end(sim)`.

## fireSense interoperability

* Accept differently-named fireSense objects: `flammableRTM` as an alternative to
  `flammableMap`, and `nonForest_timeSinceDisturbance` as an alternative to `rstTimeSinceFire`,
  with matching fallback logic.
* Declare `loadOrder = list(after = c("fireSense", "LandMine", "scfmSpread"))` so the module
  runs after any of the fire modules; moved input-object resolution from `.inputObjects` into
  `Init`.

## Fire-history retrieval and outputs

* Introduced a historical cumulative (mean-annual) burn map derived from NFDB fire polygons.
* Replaced the failing `prepInputs()` NFDB fetch with a manual `download.file()` +
  `archive::archive_extract()`, dropping invalid geometries via `terra::is.valid` (avoiding
  slow `makeValid`) and combining with `tidyterra::bind_spat_rows`; added `archive` and
  `purrr` to `reqdPkgs`.

## Housekeeping

* Added explicit `ggplot2::` and `data.table::` package prefixes throughout.
* Replaced the hand-rolled "undefined event type" warning with `noEventWarning(sim)`.
* Added LandMine, scfm, and fireSense citations.

# burnSummaries 0.0.1 (2023-09-04)

* initial module version
