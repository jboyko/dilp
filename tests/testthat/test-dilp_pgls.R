# Helper: McAbeeExample with taxonomy columns added
make_pgls_data <- function() {
  set.seed(42)
  d         <- McAbeeExample
  n         <- nrow(d)
  sp        <- paste0("Genus_sp", sample(1:5, n, replace = TRUE))
  d$species <- sp
  d$genus   <- sub("_.*", "", sp)
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
                         "placement_log", "processed_leaf_data"),
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
  expect_true(all(c("site", "species", "age_ma", "MAT.PIP", "MAP.PIP") %in% names(sp_pred)))
  expect_equal(nrow(sp_pred), 10L)  # 5 species × 2 sites
})

test_that("placement_log records every unique species", {
  result <- dilp_pgls(make_pgls_data())
  log    <- result$placement_log
  expect_s3_class(log, "data.frame")
  expect_true(all(c("species", "placement_level", "target") %in% names(log)))
  expect_equal(nrow(log), 5L)  # 5 unique species
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
  expect_equal(result$results$MAT.PIP[1], 11.25678, tolerance = 1e-4)
  expect_equal(result$results$MAP.PIP[1], 134.8828, tolerance = 1e-2)
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
  d         <- make_pgls_data()
  d$species <- NULL
  expect_error(dilp_pgls(d), "missing required columns.*species")
})

test_that("missing multiple taxonomy columns are all reported in the error", {
  d        <- make_pgls_data()
  d$genus  <- NULL
  d$order  <- NULL
  expect_error(dilp_pgls(d), "genus")
})
