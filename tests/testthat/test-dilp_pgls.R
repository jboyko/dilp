# Helper: McAbeeExample with taxonomy columns added. A species epithet belongs
# to a morphotype, not to an individual specimen, so it is assigned per
# site-morphotype; five epithets are reused to give 5 species x 2 sites.
make_pgls_data <- function() {
  d  <- McAbeeExample
  mt <- unique(paste(d$Site, d$Morphotype))
  epithets <- stats::setNames(paste0("sp", (seq_along(mt) %% 5) + 1), mt)
  d$species <- unname(epithets[paste(d$Site, d$Morphotype)])
  d$genus   <- "Betula"
  d$family  <- "Betulaceae"
  d$order   <- "Fagales"
  d$age_ma  <- 50.0
  d
}

# ----------------------------------------------------------------------
# Output structure
# ----------------------------------------------------------------------

test_that("dilp_pgls returns a list with the correct elements", {
  result <- dilp_pgls(make_pgls_data())
  expect_type(result, "list")
  expect_named(result, c("results", "species_predictions",
                         "placement_log", "unmatched_taxa",
                         "imputation_summary", "processed_leaf_data",
                         "processed_site_data"),
               ignore.order = TRUE)
})

test_that("results has one row per site with correct columns", {
  result  <- dilp_pgls(make_pgls_data())
  results <- result$results
  expect_s3_class(results, "data.frame")
  expect_true(all(c("site", "MAT.PIP", "MAP.PIP", "n_species") %in% names(results)))
  expect_equal(nrow(results), 2L)  # McAbeeExample has 2 sites
})

test_that("species_predictions has one row per species-site pair", {
  result <- dilp_pgls(make_pgls_data())
  sp_pred <- result$species_predictions
  expect_true(all(c("fossil_name", "site", "species", "age_ma",
                    "MAT.PIP", "MAP.PIP") %in% names(sp_pred)))
  expect_equal(nrow(sp_pred), 10L)  # 5 species × 2 sites
})

test_that("placement_log records every species-site occurrence", {
  result <- dilp_pgls(make_pgls_data())
  log    <- result$placement_log
  expect_s3_class(log, "data.frame")
  expect_true(all(c("fossil_name", "species", "site", "age_ma",
                    "placement_level", "target") %in% names(log)))
  expect_equal(nrow(log), 10L)  # 5 species × 2 sites
})

test_that("repeated species use each site's age and a distinct tip", {
  d <- make_pgls_data()
  sites <- unique(d$Site)
  d$age_ma <- ifelse(d$Site == sites[1], 40, 60)
  result <- dilp_pgls(d)

  expect_equal(nrow(result$placement_log), 10L)
  expect_equal(
    sort(unique(result$placement_log$age_ma)),
    c(40, 60)
  )
  expect_equal(length(unique(result$placement_log$fossil_name)), 10L)

  repeated <- result$species_predictions[
    result$species_predictions$species == result$species_predictions$species[1],
  ]
  expect_equal(nrow(repeated), 2L)
  expect_equal(sort(repeated$age_ma), c(40, 60))
  expect_equal(length(unique(repeated$fossil_name)), 2L)
})

# ----------------------------------------------------------------------
# Prediction values
# ----------------------------------------------------------------------

test_that("MAP.PIP is positive (correctly back-transformed from log scale)", {
  result <- dilp_pgls(make_pgls_data())
  expect_true(all(result$results$MAP.PIP > 0))
})

test_that("MAT.PIP is within a plausible range", {
  result <- dilp_pgls(make_pgls_data())
  expect_true(all(result$results$MAT.PIP > -10 & result$results$MAT.PIP < 50))
})

test_that("predictions are numerically stable (regression test)", {
  result <- dilp_pgls(make_pgls_data())
  expect_equal(result$results$MAT.PIP[1], 12.00846, tolerance = 1e-4)
  expect_equal(result$results$MAP.PIP[1], 118.5597, tolerance = 1e-2)
})

test_that("predictions do not depend on fossil input order", {
  d <- make_pgls_data()
  result_forward <- dilp_pgls(d)$species_predictions
  set.seed(99)
  result_shuffled <- dilp_pgls(d[sample(seq_len(nrow(d))), ])$species_predictions

  sort_predictions <- function(x) {
    x[order(x$site, x$species),
      c("site", "species", "age_ma", "MAT.PIP", "MAP.PIP")]
  }
  expect_equal(
    unname(sort_predictions(result_forward)),
    unname(sort_predictions(result_shuffled)),
    tolerance = 1e-10
  )
})

# ----------------------------------------------------------------------
# Tooth trait filling
# ----------------------------------------------------------------------

test_that("entirely untoothed assemblage produces non-NA predictions", {
  d              <- make_pgls_data()
  d[["Margin"]]  <- 1
  d[["no. primary teeth"]] <- NA
  result <- suppressWarnings(dilp_pgls(d))
  expect_false(any(is.na(result$results$MAT.PIP)))
  expect_false(any(is.na(result$results$MAP.PIP)))
})

# ----------------------------------------------------------------------
# Phylogenetic placement fallbacks
# ----------------------------------------------------------------------

test_that("family-level placement is used when genus is unknown", {
  d         <- make_pgls_data()
  d$genus   <- "UnknownGenus"
  result    <- suppressWarnings(dilp_pgls(d))
  expect_true(all(result$placement_log$placement_level == "family"))
})

test_that("order-level placement is used when genus and family are unknown", {
  d         <- make_pgls_data()
  d$genus   <- "UnknownGenus"
  d$family  <- "UnknownFamily"
  result    <- suppressWarnings(dilp_pgls(d))
  expect_true(all(result$placement_log$placement_level == "order"))
})

# ----------------------------------------------------------------------
# Input validation
# ----------------------------------------------------------------------

test_that("missing taxonomy column raises an informative error", {
  d       <- make_pgls_data()
  d$genus <- NULL
  expect_error(dilp_pgls(d), "missing required columns.*genus")
})

test_that("missing age_ma raises an informative error", {
  d        <- make_pgls_data()
  d$age_ma <- NULL
  expect_error(dilp_pgls(d), "missing required columns.*age_ma")
})

test_that("non-finite age_ma raises an informative error", {
  d        <- make_pgls_data()
  d$age_ma <- NA_real_
  expect_error(dilp_pgls(d), "age_ma must be")
})

test_that("missing multiple taxonomy columns are all reported in the error", {
  d        <- make_pgls_data()
  d$genus  <- NULL
  d$order  <- NULL
  expect_error(dilp_pgls(d), "genus")
})

# ----------------------------------------------------------------------
# Regression tests for fixed defects
# ----------------------------------------------------------------------

test_that("occurrences with no usable taxonomy are placed at the root, not fatal", {
  # An unplaceable occurrence used to be grafted under an empty tip label,
  # which corrupted the genus lookup and killed the next occurrence with a
  # zero-length node ("argument is of length zero").
  d <- make_pgls_data()
  blank <- d$Morphotype %in% head(unique(d$Morphotype), 3)
  d$species[blank] <- ""
  d$genus[blank]   <- ""
  d$family[blank]  <- ""
  d$order[blank]   <- ""

  result <- suppressWarnings(dilp_pgls(d))
  expect_true(all(result$placement_log$placed))
  expect_true(any(result$placement_log$placement_level == "root"))
  expect_true(all(is.finite(result$results$MAT.PIP)))
  expect_true(all(is.finite(result$results$MAP.PIP)))
})

test_that("unmatched taxa are reported for the caller to clean up", {
  d        <- make_pgls_data()
  d$genus  <- "Mispelling"
  d$family <- "Mispellingaceae"
  d$order  <- "Mispellingales"

  result <- suppressWarnings(dilp_pgls(d))
  expect_s3_class(result$unmatched_taxa, "data.frame")
  expect_equal(nrow(result$unmatched_taxa), nrow(result$species_predictions))
  expect_true(all(result$unmatched_taxa$genus_reported == "Mispelling"))
})

test_that("missing predictors are imputed inside the function", {
  d <- make_pgls_data()
  d[["Petiole width"]] <- NA

  result <- suppressWarnings(dilp_pgls(d))
  summary <- result$imputation_summary
  expect_s3_class(summary, "data.frame")
  expect_gt(summary$n_imputed[summary$predictor == "petiole_metric"], 0)
  expect_true(all(is.finite(result$species_predictions$MAT.PIP)))
  expect_true(all(is.finite(result$species_predictions$MAP.PIP)))
})

test_that("unnamed morphotypes stay distinct occurrences", {
  # Grouping by taxonomy alone collapsed every unnamed morphotype at a site
  # into a single row, silently merging distinct fossils.
  d          <- make_pgls_data()
  d$species  <- NA_character_
  n_expected <- nrow(unique(d[, c("Site", "Morphotype")]))

  result <- suppressWarnings(dilp_pgls(d))
  expect_equal(nrow(result$species_predictions), n_expected)
  expect_setequal(result$species_predictions$species, unique(d$Morphotype))
})

test_that("site MAP is the geometric mean of occurrence predictions", {
  result <- dilp_pgls(make_pgls_data())
  sp     <- result$species_predictions
  site   <- result$results$site[1]
  expect_equal(
    result$results$MAP.PIP[1],
    exp(mean(log(sp$MAP.PIP[sp$site == site]))),
    tolerance = 1e-10
  )
})

test_that("informal (quoted) taxonomy is censored unless asked for", {
  d       <- make_pgls_data()
  d$genus <- '"Betula"'

  formal   <- suppressWarnings(dilp_pgls(d, taxonomy_scenario = "formal_only"))
  informal <- suppressWarnings(dilp_pgls(d, taxonomy_scenario = "include_informal"))

  expect_true(all(formal$placement_log$placement_level == "family"))
  expect_true(all(informal$placement_log$placement_level == "genus"))
})
