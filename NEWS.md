# burnSummaries (development version)

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
