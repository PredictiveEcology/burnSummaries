# burnSummaries (development version)

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
