# dilp 1.4.0.9001
* `dilp_pgls()` now takes raw leaf-level data all the way through, so no
  pre-processing is required of the caller. Specimens are averaged to
  morphotypes, morphotypes to species-by-site occurrences, and missing
  predictors are imputed inside the function with the bagged-tree imputers
  fitted on the calibration data. Values imputed by the caller beforehand were
  previously discarded by `dilp_processing()` without warning.
* Fixed a fatal error ("argument is of length zero") when an occurrence had no
  usable taxonomy. Such an occurrence was grafted under an empty tip label,
  which put an `NA` into the genus lookup and made the *next* occurrence
  resolve to a zero-length node. Taxonomic lookups are now restricted to an
  immutable snapshot of the scaffold tips, so a grafted fossil can never
  become evidence for a later placement, and a failed lookup falls through to
  the next rank instead of producing a zero-length index.
* Fixed occurrence keying: occurrences were grouped by taxonomy alone, so every
  morphotype without a binomial at a site collapsed into a single row. An
  occurrence is now identified by its binomial where genus and epithet are
  both given and by its morphotype otherwise.
* Fixed the phylogenetic correction, which omitted the fitted lambda scaling of
  the fossil-to-calibration covariance, and which pruned the tree with
  `keep.tip()` before computing covariances — pruning can move the root and
  rescale every entry. The covariance matrix of the grafted tree is now
  subsetted instead.
* Site-level MAP is now a geometric mean, averaged on the log scale the model
  was fitted on and then back-transformed, rather than an arithmetic mean of
  back-transformed values.
* `age_ma` is validated, and is required for every row.
* Added the `taxonomy_scenario` argument. Informal taxonomic names, marked with
  quotation marks, are censored under the default `"formal_only"` and used as
  provisional placement hypotheses under `"include_informal"`.
* Added `$unmatched_taxa`, listing occurrences whose reported taxonomy matched
  nothing in the calibration phylogeny below the root, to help find
  misspellings and synonyms; and `$imputation_summary`, giving the number of
  imputed values per predictor.
* Refreshed the pre-fitted PGLS model objects, scaffold tree, and RMSE
  constants from the current calibration pipeline. `dilp_pgls()` now reproduces
  the analysis pipeline's fossil predictions to within 5e-13 degrees C and
  2e-11 cm across all 360 occurrences in both taxonomy scenarios.
* Added `caret` and `ipred` as package dependencies, for imputation.

# dilp 1.4.0
* Added `dilp_pgls()` for Phylogenetically-Informed Predictions (PIP). Accepts
  the same specimen-level input as `dilp()` plus taxonomy columns (`species`,
  `genus`, `family`, `order`) and fossil age (`age_ma`). Returns site-level MAT
  and MAP estimates with a phylogenetic covariance correction, per-species
  predictions, and a placement log recording how each fossil taxon was grafted
  onto the angiosperm scaffold tree.
* Added `ape` and `phytools` as package dependencies.

# dilp 1.3.0
* Refixed outliers merge to improve readability
* Added some extra visiblity parameters in `dilp_cca()` that still don't completely work

# dilp 1.2.1
* Fixed Outliers merge so there aren't thousands of outliers.
* Adjusted Errors so that it pulls by column header rather than index

# dilp 1.2.0
* Added Miller(2006) error calculations to `temp_slr()`
* Added KowalskiDilcher (2003) regression to `temp_slr()`
* Added `dilp_whittaker()` function for plotting localities to biome.
* Makeover for error and outlier reporting from `dilp()`

# dilp 1.1.0
* Added `view_regressions()` to allow users to easily see their options for preloaded regressions in dilp, lma, temp, and precip.
* Fixed `dilp()` so that the mixed_margins aren't screwing everything up.
* Fixed `dilp()` documentation to properly reference the correct parameter options.
* Entry field "no_primary_teeth" changed to "no_of_primary_teeth" for consistency
* Updated provided data entry template so that column names match listed expected column names
* Added conditionals to `dilp()` so that if you enter blanks for raw blade area, and raw blade perimeter, but length of cut perimeter is 0, blade perimeter/blade area are present, and the leaf is toothed, the raw values are automatically filled in with the blade perimeter/blade area values.
* Added conditionals to `dilp()` site_mean values so that a locality with only untoothed leaves is assigned a site tc_ip of 0 and a site perimeter_ratio of 1, allowing MAT estimates to be formed.

# dilp 1.0.0
* Initial CRAN submission
* `dilp()` given new argument `subsite_cols` to allow multiple tiers of specimen/morphotype organization.
* `dilp()` given new options for regression parameters in `params`.  Choose between "global" and "northern_hemisphere".

# dilp (development version 0.0.0.9001)
* Functions implemented for simple linear regression (temp_slr() and precip_slr())

# dilp (development version 0.0.0.9000)

* Functions implemented for digital leaf physiognomy (dilp())
* Functions implemented for leaf mass per area reconstruction (lma())
