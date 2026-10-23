# Trait columns (dilp-style names) used in the PGLS model.
# Must match the columns computed by dilp_processing().
.pgls_trait_cols <- c(
  "petiole_metric", "ln_leaf_area", "fdr", "margin",
  "perimeter_ratio", "tc_p", "tc_ip", "avg_ta",
  "ta_ba", "ta_p", "ta_ip", "tc_ba"
)

# Tooth-derived traits are structural zeros on an untoothed leaf.
.pgls_tooth_zero_cols <- c("tc_p", "tc_ip", "avg_ta", "ta_ba", "ta_p", "ta_ip", "tc_ba")

# The shipped imputers were fitted on the calibration data under its own
# trait names, so fossil columns are renamed before and after imputation.
.pgls_dilp_to_calibration <- c(
  petiole_metric  = "pw2.a.ratio",
  ln_leaf_area    = "ln.leaf.area.mm2",
  fdr             = "feret.diam.ratio",
  margin          = "margin.score",
  perimeter_ratio = "perim.ratio",
  tc_p            = "teeth.perimeter.percm",
  tc_ip           = "teeth.interior.percm",
  avg_ta          = "avt.tooth.area",
  ta_ba           = "tooth.area.blade.area.ratio",
  ta_p            = "tooth.area.perimeter",
  ta_ip           = "tooth.area.interior",
  tc_ba           = "teeth.blade.area.ratio"
)

#' Generate phylogenetically-informed DiLP climate reconstructions
#'
#' @description
#' `dilp_pgls()` uses Phylogenetically-Informed Predictions (PIP) to estimate
#' mean annual temperature (MAT) and mean annual precipitation (MAP) from fossil
#' leaf specimens. Predictions borrow statistical strength from phylogenetically
#' close extant calibration species, improving estimates particularly when few
#' specimens are available per site.
#'
#' The method fits a Phylogenetic Generalised Least Squares (PGLS) model on a
#' large extant training dataset (Royer et al., n = 1847 species), then at
#' prediction time adds a phylogenetic covariance correction that pulls each
#' fossil occurrence's prediction toward the climate values of its closest
#' extant relatives on the tree.
#'
#' Leaf-level input is taken all the way through: specimens are averaged to
#' morphotypes, morphotypes to species-by-site occurrences, missing predictors
#' are imputed with the bagged-tree imputers fitted on the calibration data,
#' and each occurrence is grafted onto the scaffold phylogeny at its site's
#' age. No pre-processing or imputation is required of the caller.
#'
#' @param specimen_data A data frame of specimen-level (leaf-level) leaf
#'   physiognomic data. Must include all columns required by
#'   [dilp_processing()], plus:
#'   * `genus`, `family`, `order` - taxonomic ranks used for phylogenetic
#'     placement, in that order of preference
#'   * `age_ma` - site age in millions of years, used for occurrence-tip depth
#'
#'   An optional `species` column holds the species epithet. Occurrences are
#'   identified by `genus` plus `species` where both are given, and by
#'   `morphotype` otherwise, so unnamed morphotypes stay distinct.
#'
#'   Placement proceeds genus -> family -> order -> root. Occurrences that
#'   reach only the root, and any that fail placement outright, are reported in
#'   `$unmatched_taxa`.
#'
#' @param taxonomy_scenario How to treat informal taxonomic names, which Dana
#'   Royer's datasets mark with quotation marks (for example `"Cornophyllum"`).
#'   `"formal_only"` (the default) censors every quoted rank, so an informal
#'   name is never used as phylogenetic evidence. `"include_informal"` strips
#'   the quotation marks and uses the reported name as a provisional placement
#'   hypothesis; it is a sensitivity analysis, not a second best estimate.
#'
#' @return A list with seven elements:
#' * `results` - a data frame of site-level MAT and MAP estimates with columns
#'   `site`, `MAT.PIP` (degrees C), `MAP.PIP` (cm), `n_species`, `MAT.PIP.error`
#'   (+/- degrees C), `MAP.PIP.error.plus`, and `MAP.PIP.error.minus` (cm). MAP is a
#'   geometric site mean, averaged on the log scale and back-transformed.
#'   Uncertainty columns are RMSE-based, from 10-fold site-grouped
#'   cross-validation in Boyko et al. (in prep), and do not reflect
#'   phylogenetic placement uncertainty.
#' * `species_predictions` - a data frame of per-occurrence estimates with
#'   columns `fossil_name`, `site`, `species`, `morphotype`, `age_ma`,
#'   `MAT.PIP`, and `MAP.PIP`.
#' * `placement_log` - a data frame recording how each occurrence was placed
#'   on the phylogeny (`fossil_name`, `species`, `site`, `age_ma`, `placed`,
#'   `placement_level`, `placement_target`, `age_fallback`, `error`).
#' * `unmatched_taxa` - a data frame of occurrences whose reported taxonomy
#'   matched nothing in the calibration phylogeny, with the names as reported
#'   (`genus_reported`, `family_reported`, `order_reported`) and whether each
#'   rank was informal. Use it to find misspellings and synonyms.
#' * `imputation_summary` - a data frame with one row per predictor giving the
#'   number of occurrences whose value was imputed.
#' * `processed_leaf_data` - the full processed specimen-level data frame from
#'   [dilp_processing()].
#' * `processed_site_data` - site-level means of leaf physiognomic variables,
#'   aggregated via morphotype then site (same structure as [dilp()]). Used by
#'   [dilp_cca()].
#'
#' @references
#' * Boyko, J.D., Hagen, E.R., O'Meara, B.C., Lawing, A.M., Peppe, D.J.,
#'   and Royer, D.L. (in prep). Phylogenetically-Informed DiLP:
#'   paleoclimate reconstruction from fossil leaves using PGLS and
#'   Phylogenetically-Informed Predictions.
#' * Freckleton, R.P. (2015). PIC computation in O(N) time.
#'   In Garamszegi, L.Z. (Ed.), \emph{Modern Phylogenetic Comparative Methods
#'   and Their Application in Evolutionary Biology}. Springer.
#' * Gardner, J.D., Baker, J., Venditti, C., and Organ, C.L. (2024).
#'   Phylogenetically-Informed Predictions. \emph{Methods in Ecology
#'   and Evolution}.
#' * Lowe, A.J., Flynn, A.G., Butrim, M.J., Baumgartner, A., Peppe, D.J.,
#'   and Royer, D.L. (2024). Reconstructing terrestrial paleoclimate and
#'   paleoecology with fossil leaves using Digital Leaf Physiognomy and leaf
#'   mass per area. \emph{J. Vis. Exp.} (212), e66838.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Add required taxonomy and age columns to McAbeeExample
#' fossil_data <- McAbeeExample
#' fossil_data$genus  <- "Betula"
#' fossil_data$family <- "Betulaceae"
#' fossil_data$order  <- "Fagales"
#' fossil_data$age_ma <- 52.0
#'
#' results <- dilp_pgls(fossil_data)
#' results$results
#' results$species_predictions
#' results$unmatched_taxa
#'
#' # Plot sites on a Whittaker biome diagram
#' dilp_whittaker(results)
#'
#' # Check fossil site physiognomy against the calibration space
#' dilp_cca(results)
#' }
dilp_pgls <- function(specimen_data,
                      taxonomy_scenario = c("formal_only",
                                            "include_informal")) {

  taxonomy_scenario <- match.arg(taxonomy_scenario)

  # ------------------------------------------------------------------
  # 1. Standardise column names and check what the caller supplied
  # ------------------------------------------------------------------
  colnames(specimen_data) <- colnameClean(specimen_data)

  required_cols <- c("site", "morphotype", "genus", "family", "order", "age_ma")
  missing_cols  <- required_cols[!required_cols %in% colnames(specimen_data)]
  if (length(missing_cols) > 0) {
    stop(paste("specimen_data is missing required columns:",
               paste(missing_cols, collapse = ", ")))
  }
  if (!"species" %in% colnames(specimen_data)) {
    specimen_data$species <- NA_character_
  }

  specimen_data$site       <- trimws(as.character(specimen_data$site))
  specimen_data$morphotype <- trimws(as.character(specimen_data$morphotype))
  if (anyNA(specimen_data$site) || any(!nzchar(specimen_data$site))) {
    stop("Every row of specimen_data needs a site.")
  }
  if (anyNA(specimen_data$morphotype) || any(!nzchar(specimen_data$morphotype))) {
    stop("Every row of specimen_data needs a morphotype.")
  }

  specimen_data$age_ma <- suppressWarnings(as.numeric(specimen_data$age_ma))
  if (any(!is.finite(specimen_data$age_ma)) || any(specimen_data$age_ma < 0)) {
    stop("age_ma must be a finite, non-negative site age for every row.")
  }
  age_spread <- tapply(specimen_data$age_ma, specimen_data$site,
                       function(x) diff(range(x)))
  if (any(age_spread > 0)) {
    warning("age_ma varies within site(s): ",
            paste(names(age_spread)[age_spread > 0], collapse = ", "),
            ". Occurrence ages are averaged within site-morphotype.",
            call. = FALSE)
  }

  # ------------------------------------------------------------------
  # 2. Process leaf traits
  # ------------------------------------------------------------------
  processed <- dilp_processing(specimen_data)

  # dilp_processing() drops no columns, so taxonomy survives

  # ------------------------------------------------------------------
  # 3. Average specimens to morphotypes exactly as dilp() does, so the
  #    mixed-margin rule and tooth conventions match the calibration data.
  # ------------------------------------------------------------------
  site_morphotype <- processed %>%
    dplyr::select(-"specimen_number") %>%
    dplyr::summarize(
      .by = c("site", "morphotype"),
      dplyr::across(dplyr::where(~ !is.character(.)), \(x) mean(x, na.rm = TRUE))
    )

  # Mixed-margin morphotypes get margin = 0.5
  if (length(unique(stats::na.omit(site_morphotype$margin))) > 2) {
    site_morphotype$margin[site_morphotype$margin > 0 & site_morphotype$margin < 1] <- 0.5
  }

  processed_site_data <- site_morphotype %>%
    dplyr::group_by(.data$site) %>%
    dplyr::select(-"morphotype") %>%
    dplyr::summarize(dplyr::across(dplyr::where(~ !is.character(.)), \(x) mean(x, na.rm = TRUE))) %>%
    dplyr::mutate(
      margin         = .data$margin * 100,
      tc_ip          = ifelse(.data$margin == 100, 0,  .data$tc_ip),
      perimeter_ratio = ifelse(.data$margin == 100, 1, .data$perimeter_ratio)
    )

  # Averaging an all-untoothed morphotype gives NaN for every tooth trait.
  # Those are biologically real zeros, not missing data, and must be set
  # before imputation so the imputer never invents a tooth measurement.
  site_morphotype <- as.data.frame(site_morphotype)
  untoothed <- !is.na(site_morphotype$margin) & site_morphotype$margin == 1
  for (col in .pgls_tooth_zero_cols) {
    site_morphotype[[col]][untoothed & !is.finite(site_morphotype[[col]])] <- 0
  }
  site_morphotype$perimeter_ratio[
    untoothed & !is.finite(site_morphotype$perimeter_ratio)
  ] <- 1

  # ------------------------------------------------------------------
  # 4. Attach the taxonomy reported for each site-morphotype
  # ------------------------------------------------------------------
  taxonomy <- .pgls_morphotype_taxonomy(specimen_data)
  taxonomy <- taxonomy[match(
    paste(site_morphotype$site, site_morphotype$morphotype, sep = "\r"),
    paste(taxonomy$site, taxonomy$morphotype, sep = "\r")
  ), ]
  if (anyNA(taxonomy$site)) {
    stop("Failed to reattach taxonomy to one or more processed morphotypes.")
  }

  # An occurrence is identified by its binomial where one was assigned, and by
  # its morphotype otherwise. Collapsing every unnamed morphotype at a site
  # into a single "species" would silently merge distinct fossils.
  genus_label     <- .pgls_strip_quotes(taxonomy$genus)
  species_epithet <- .pgls_strip_quotes(taxonomy$species)
  genus_label[is.na(genus_label)]         <- ""
  species_epithet[is.na(species_epithet)] <- ""
  named_taxon <- gsub("[[:space:]]+", " ",
                      trimws(paste(genus_label, species_epithet)))
  site_morphotype$species <- ifelse(
    nzchar(genus_label) & nzchar(species_epithet),
    named_taxon,
    site_morphotype$morphotype
  )
  site_morphotype$genus_reported  <- .pgls_normalise_taxon(taxonomy$genus)
  site_morphotype$family_reported <- .pgls_normalise_taxon(taxonomy$family)
  site_morphotype$order_reported  <- .pgls_normalise_taxon(taxonomy$order)

  # ------------------------------------------------------------------
  # 5. Aggregate morphotypes within species-by-site occurrences. Traits and
  #    ages are never pooled across sites.
  # ------------------------------------------------------------------
  mean_or_na <- function(x) {
    value <- mean(x, na.rm = TRUE)
    if (is.nan(value)) NA_real_ else value
  }

  keys <- site_morphotype[, c("species", "site")]
  sp_site <- stats::aggregate(
    site_morphotype[, c("age_ma", .pgls_trait_cols)],
    by = keys, FUN = mean_or_na
  )
  morphotypes <- stats::aggregate(
    list(morphotype = site_morphotype$morphotype), by = keys,
    FUN = function(x) paste(sort(unique(x)), collapse = "|")
  )
  reported <- stats::aggregate(
    site_morphotype[, c("genus_reported", "family_reported", "order_reported")],
    by = keys, FUN = .pgls_resolve_reported_taxon
  )
  sp_site <- merge(sp_site, morphotypes, by = c("species", "site"), sort = FALSE)
  sp_site <- merge(sp_site, reported, by = c("species", "site"), sort = FALSE)

  sp_site <- .pgls_apply_taxonomy_scenario(sp_site, taxonomy_scenario)
  sp_site$fossil_name <- make.unique(gsub(
    "[^[:alnum:]]+", "_", paste(sp_site$species, sp_site$site)
  ))

  # ------------------------------------------------------------------
  # 6. Impute missing predictors with the calibration-fitted imputers.
  #    Imputation happens here, inside the function, so that callers never
  #    have to do it themselves -- and so that values they impute upstream
  #    cannot be silently discarded by dilp_processing() above.
  # ------------------------------------------------------------------
  imputed <- .pgls_impute_traits(sp_site)
  imputation_summary <- imputed$summary
  if (any(imputation_summary$n_imputed > 0)) {
    filled <- imputation_summary[imputation_summary$n_imputed > 0, ]
    message(
      "Imputed missing predictors for ", sum(filled$n_imputed),
      " occurrence-trait values using the calibration imputers:\n  ",
      paste(paste0(filled$predictor, " (", filled$n_imputed, ")"),
            collapse = ", ")
    )
  }

  # ------------------------------------------------------------------
  # 7. Load scaffold tree and place fossil occurrences
  # ------------------------------------------------------------------
  scaffold_path <- system.file("extdata", "tre_scaffold.tre", package = "dilp")
  if (!nzchar(scaffold_path)) {
    stop("Scaffold tree not found. Please reinstall the dilp package.")
  }
  tree <- ape::read.tree(scaffold_path)

  placement     <- .pgls_place_fossils(tree, sp_site, pgls_name_table)
  tree          <- placement$tree
  placement_log <- placement$log

  unmatched <- placement_log$placed & placement_log$placement_level == "root"
  if (any(unmatched)) {
    warning(sum(unmatched), " occurrence(s) matched no calibration taxon ",
            "below the root and were placed at the root. See $unmatched_taxa.",
            call. = FALSE)
  }
  if (any(!placement_log$placed)) {
    failed <- placement_log[!placement_log$placed, ]
    warning("Dropping ", nrow(failed), " occurrence(s) that failed ",
            "phylogenetic placement: ",
            paste(failed$fossil_name, collapse = ", "), call. = FALSE)
  }

  placed_names <- placement_log$fossil_name[placement_log$placed]
  keep         <- sp_site$fossil_name %in% placed_names
  unmatched_taxa <- .pgls_unmatched_taxa(sp_site, placement_log)
  sp_site      <- sp_site[keep, , drop = FALSE]
  trait_mat    <- imputed$mat[keep, , drop = FALSE]
  trait_map    <- imputed$map[keep, , drop = FALSE]

  if (nrow(sp_site) == 0) {
    stop("No fossil occurrences could be placed on the phylogeny.")
  }

  # ------------------------------------------------------------------
  # 8. Phylogenetic adjustment: lambda * t(V_cross) %*% Ke
  #    Subset the covariance matrix of the grafted tree rather than pruning
  #    the tree: dropping tips can move the root and rescale the covariance.
  # ------------------------------------------------------------------
  fossil_occurrences <- sp_site$fossil_name
  idx_extant <- union(names(pgls_Ke_mat), names(pgls_Ke_map))
  phylomat <- ape::vcv(tree)
  phylomat <- phylomat[c(idx_extant, fossil_occurrences),
                       c(idx_extant, fossil_occurrences), drop = FALSE]

  phylo_adj <- function(Ke, lambda) {
    V_cross <- phylomat[names(Ke), fossil_occurrences, drop = FALSE]
    stats::setNames(lambda * as.numeric(t(V_cross) %*% Ke), fossil_occurrences)
  }
  phylo_adj_mat <- phylo_adj(pgls_Ke_mat, pgls_lambda_mat)
  phylo_adj_map <- phylo_adj(pgls_Ke_map, pgls_lambda_map)

  # ------------------------------------------------------------------
  # 9. Build design matrices and predict per occurrence
  # ------------------------------------------------------------------
  design <- function(traits) {
    X <- cbind(`(Intercept)` = 1, as.matrix(traits))
    X[, names(pgls_beta_mat), drop = FALSE]
  }

  sp_site$MAT.PIP <- as.numeric(design(trait_mat) %*% pgls_beta_mat) +
                       phylo_adj_mat[sp_site$fossil_name]
  log_map_pip     <- as.numeric(design(trait_map) %*% pgls_beta_map) +
                       phylo_adj_map[sp_site$fossil_name]
  sp_site$MAP.PIP <- exp(log_map_pip)

  # ------------------------------------------------------------------
  # 10. Aggregate occurrence predictions to site level. MAT averages
  #     arithmetically; MAP averages on the log scale it was fitted on and is
  #     then back-transformed, giving a geometric site mean.
  # ------------------------------------------------------------------
  sp_site$.log_map <- log_map_pip
  results <- sp_site %>%
    dplyr::group_by(.data$site) %>%
    dplyr::summarize(
      MAT.PIP   = mean(.data$MAT.PIP, na.rm = TRUE),
      MAP.PIP   = exp(mean(.data$.log_map, na.rm = TRUE)),
      n_species = dplyr::n(),
      .groups   = "drop"
    )
  sp_site$.log_map <- NULL

  # ------------------------------------------------------------------
  # 11. Attach RMSE-based uncertainty from 10-fold site-grouped
  #     cross-validation (Boyko et al., in prep). MAT error is symmetric
  #     (degrees C); MAP errors are asymmetric, back-transformed from the log scale.
  # ------------------------------------------------------------------
  log_map_site <- log(results$MAP.PIP)
  results$MAT.PIP.error       <- pgls_rmse_mat
  results$MAP.PIP.error.plus  <- exp(log_map_site + pgls_rmse_log_map) - results$MAP.PIP
  results$MAP.PIP.error.minus <- results$MAP.PIP - exp(log_map_site - pgls_rmse_log_map)

  message(
    "Uncertainty columns (MAT.PIP.error, MAP.PIP.error.*) are RMSE-based, ",
    "derived from 10-fold site-grouped cross-validation reported in ",
    "Boyko et al. (in prep). They do not account for phylogenetic placement uncertainty."
  )

  list(
    results             = as.data.frame(results),
    species_predictions = as.data.frame(
      sp_site[, c("fossil_name", "site", "species", "morphotype", "age_ma",
                  "MAT.PIP", "MAP.PIP")]
    ),
    placement_log       = placement_log,
    unmatched_taxa      = unmatched_taxa,
    imputation_summary  = imputation_summary,
    processed_leaf_data = processed,
    processed_site_data = as.data.frame(processed_site_data)
  )
}


# ----------------------------------------------------------------------
# Internal: informal (quoted) taxonomy handling
#
# Dana Royer (pers. comm., April 2026) specified that order, family, or genus
# names inside quotation marks are informal and must not be used as
# phylogenetic evidence in the primary analysis.
# ----------------------------------------------------------------------
# Straight and curly quotation marks, built from code points to keep this
# file ASCII-only: " and the four Unicode marks U+2018, U+2019, U+201C, U+201D.
.pgls_quote_marks <- c('"', vapply(c(8216L, 8217L, 8220L, 8221L),
                                   intToUtf8, character(1)))

.pgls_is_informal <- function(x) {
  x <- as.character(x)
  informal <- rep(FALSE, length(x))
  for (mark in .pgls_quote_marks) {
    informal <- informal | grepl(mark, x, fixed = TRUE)
  }
  informal & !is.na(x)
}

.pgls_strip_quotes <- function(x) {
  x <- as.character(x)
  for (mark in .pgls_quote_marks) {
    x <- gsub(mark, "", x, fixed = TRUE)
  }
  trimws(gsub("[[:space:]]+", " ", x))
}

.pgls_normalise_taxon <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | !nzchar(x) | tolower(x) == "unknown"] <- "unknown"
  x
}

.pgls_resolve_reported_taxon <- function(x) {
  x <- .pgls_normalise_taxon(x)
  known <- unique(x[x != "unknown"])
  if (length(known) > 1) {
    stop("Conflicting non-missing taxonomy within a species-site: ",
         paste(known, collapse = ", "))
  }
  if (length(known) == 1) known else "unknown"
}

.pgls_apply_taxonomy_scenario <- function(data, scenario) {
  for (rank in c("genus", "family", "order")) {
    reported <- .pgls_normalise_taxon(data[[paste0(rank, "_reported")]])
    informal <- .pgls_is_informal(reported)
    cleaned  <- .pgls_normalise_taxon(.pgls_strip_quotes(reported))
    if (scenario == "formal_only") cleaned[informal] <- "unknown"

    data[[paste0(rank, "_reported")]] <- reported
    data[[paste0(rank, "_informal")]] <- informal
    data[[rank]] <- cleaned
  }
  data$taxonomy_scenario <- scenario
  data
}

# One taxonomy per site-morphotype. The placement ranks must be unambiguous,
# because they are the phylogenetic evidence; a morphotype carrying two
# different genera, families or orders is a data error and stops the run. A
# conflicting species epithet only affects the occurrence's label, so the
# first non-missing value is used and the caller is warned.
.pgls_morphotype_taxonomy <- function(specimen_data) {
  cols <- c("site", "morphotype", "genus", "family", "order", "species")
  records <- unique(specimen_data[, cols])
  key <- paste(records$site, records$morphotype, sep = "\r")

  rank_key <- paste(key, records$genus, records$family, records$order, sep = "\r")
  conflicted <- unique(key[duplicated(unique(cbind(key, rank_key))[, "key"])])
  if (length(conflicted)) {
    stop("Conflicting genus/family/order within site-morphotype(s): ",
         paste(sub("\r", " / ", conflicted), collapse = "; "))
  }

  first_known <- function(x) {
    known <- x[!is.na(x) & nzchar(trimws(as.character(x)))]
    if (length(known)) known[1] else NA_character_
  }
  split_rows <- split(seq_len(nrow(records)), key)
  ambiguous <- names(split_rows)[vapply(
    split_rows,
    function(i) length(unique(first_known(records$species[i]))) > 0 &&
                length(unique(stats::na.omit(records$species[i]))) > 1,
    logical(1)
  )]
  if (length(ambiguous)) {
    warning("More than one species epithet within site-morphotype(s): ",
            paste(sub("\r", " / ", ambiguous), collapse = "; "),
            ". Using the first non-missing epithet for each.", call. = FALSE)
  }

  taxonomy <- records[vapply(split_rows, `[`, integer(1), 1L), , drop = FALSE]
  taxonomy$species <- vapply(split_rows,
                             function(i) first_known(records$species[i]),
                             character(1))
  row.names(taxonomy) <- NULL
  taxonomy
}


# ----------------------------------------------------------------------
# Internal: impute missing predictors with the calibration-fitted
# bagged-tree imputers. The MAT and MAP imputers are deliberately separate,
# matching the two calibration fits.
# ----------------------------------------------------------------------
.pgls_impute_traits <- function(sp_site) {
  if (!requireNamespace("caret", quietly = TRUE) ||
      !requireNamespace("ipred", quietly = TRUE)) {
    stop("Imputation of missing leaf traits needs the caret and ipred ",
         "packages. Install them with ",
         'install.packages(c("caret", "ipred")).')
  }

  raw <- sp_site[, .pgls_trait_cols, drop = FALSE]
  colnames(raw) <- unname(.pgls_dilp_to_calibration[colnames(raw)])

  apply_imputer <- function(imputer) {
    filled <- as.data.frame(stats::predict(imputer, newdata = raw))
    filled <- filled[, unname(.pgls_dilp_to_calibration), drop = FALSE]
    colnames(filled) <- names(.pgls_dilp_to_calibration)
    out <- filled[, .pgls_trait_cols, drop = FALSE]
    if (anyNA(out)) {
      still_missing <- names(which(colSums(is.na(out)) > 0))
      stop("Imputation left missing predictors: ",
           paste(still_missing, collapse = ", "))
    }
    out
  }

  summary <- data.frame(
    predictor = .pgls_trait_cols,
    n_imputed = as.integer(colSums(is.na(sp_site[, .pgls_trait_cols]))),
    stringsAsFactors = FALSE
  )
  row.names(summary) <- NULL

  list(
    mat     = apply_imputer(pgls_impute_model_mat_traits),
    map     = apply_imputer(pgls_impute_model_map_traits),
    summary = summary
  )
}


# ----------------------------------------------------------------------
# Internal: report occurrences whose taxonomy matched nothing below the root
# ----------------------------------------------------------------------
.pgls_unmatched_taxa <- function(sp_site, placement_log) {
  unmatched <- !placement_log$placed | placement_log$placement_level == "root"
  log_rows  <- placement_log[unmatched, , drop = FALSE]
  if (nrow(log_rows) == 0) {
    return(data.frame(
      fossil_name = character(0), site = character(0), species = character(0),
      genus_reported = character(0), family_reported = character(0),
      order_reported = character(0), genus_informal = logical(0),
      family_informal = logical(0), order_informal = logical(0),
      placed = logical(0), stringsAsFactors = FALSE
    ))
  }
  idx <- match(log_rows$fossil_name, sp_site$fossil_name)
  out <- data.frame(
    fossil_name     = log_rows$fossil_name,
    site            = log_rows$site,
    species         = log_rows$species,
    genus_reported  = sp_site$genus_reported[idx],
    family_reported = sp_site$family_reported[idx],
    order_reported  = sp_site$order_reported[idx],
    genus_informal  = sp_site$genus_informal[idx],
    family_informal = sp_site$family_informal[idx],
    order_informal  = sp_site$order_informal[idx],
    placed          = log_rows$placed,
    stringsAsFactors = FALSE
  )
  row.names(out) <- NULL
  out
}


# ----------------------------------------------------------------------
# Internal: graft fossil occurrences onto scaffold tree
#
# The scaffold tip labels are an immutable snapshot. `tree` grows as each
# occurrence is grafted, but a fossil must never become taxonomic evidence for
# a later fossil's placement -- and an empty or unparsable fossil label must
# never enter the genus lookup.
# ----------------------------------------------------------------------
.pgls_place_fossils <- function(tree, occurrences, name_table) {
  scaffold_tips   <- tree$tip.label
  scaffold_genera <- sub("_.*$", "", scaffold_tips)
  log_rows <- vector("list", nrow(occurrences))

  known <- function(value) {
    length(value) == 1L && !is.na(value) && nzchar(value) && value != "unknown"
  }

  for (j in seq_len(nrow(occurrences))) {
    result <- .pgls_graft_one(
      tree, scaffold_tips, scaffold_genera, name_table,
      tip_label = occurrences$fossil_name[j],
      age_ma    = occurrences$age_ma[j],
      genus     = occurrences$genus[j],
      family    = occurrences$family[j],
      order     = occurrences$order[j],
      known     = known
    )
    tree <- result$tree
    log_rows[[j]] <- data.frame(
      fossil_name      = occurrences$fossil_name[j],
      species          = occurrences$species[j],
      site             = occurrences$site[j],
      age_ma           = occurrences$age_ma[j],
      placed           = result$placed,
      placement_level  = result$placement_level,
      placement_target = result$placement_target,
      # Retained alias: earlier versions named this column `target`.
      target           = result$placement_target,
      age_fallback     = result$age_fallback,
      error            = result$error,
      stringsAsFactors = FALSE
    )
  }

  log_df <- do.call(rbind, log_rows)
  row.names(log_df) <- NULL
  list(tree = tree, log = log_df)
}

.pgls_graft_one <- function(tree, scaffold_tips, scaffold_genera, name_table,
                            tip_label, age_ma, genus, family, order, known) {
  failure <- function(message) {
    list(tree = tree, placed = FALSE, placement_level = NA_character_,
         placement_target = NA_character_, age_fallback = NA_character_,
         error = message)
  }

  if (length(tip_label) != 1L || is.na(tip_label) || !nzchar(tip_label)) {
    return(failure("tip_label must be one non-empty string"))
  }
  if (length(age_ma) != 1L || !is.finite(age_ma) || age_ma < 0) {
    return(failure("age_ma must be one finite, non-negative value"))
  }
  if (tip_label %in% scaffold_tips) {
    return(list(tree = tree, placed = TRUE, placement_level = "existing_tip",
                placement_target = tip_label, age_fallback = "none",
                error = NA_character_))
  }

  # Every taxonomic lookup is restricted to the original scaffold.
  anchor_tips <- intersect(scaffold_tips, tree$tip.label)
  if (!length(anchor_tips)) {
    return(failure("no scaffold tips remain in tree"))
  }
  anchor_genera <- scaffold_genera[match(anchor_tips, scaffold_tips)]

  # Returns NULL -- never a zero-length index -- when nothing matches, so the
  # next fallback runs instead of a zero-length node reaching nodeheight().
  lookup_target <- function(candidate_tips) {
    candidate_tips <- intersect(anchor_tips, candidate_tips)
    if (length(candidate_tips) >= 2L) {
      node <- ape::getMRCA(tree, candidate_tips)
      return(if (length(node) == 1L && !is.na(node)) node else NULL)
    }
    if (length(candidate_tips) == 1L) {
      node <- match(candidate_tips, tree$tip.label)
      return(if (length(node) == 1L && !is.na(node)) node else NULL)
    }
    NULL
  }

  target_node      <- NULL
  placement_level  <- "root"
  placement_target <- "root"
  age_fallback     <- "none"

  if (known(genus)) {
    target_node <- lookup_target(anchor_tips[anchor_genera == genus])
    if (!is.null(target_node)) {
      placement_level  <- "genus"
      placement_target <- genus
    }
  }
  if (is.null(target_node) && known(family)) {
    family_genera <- unique(name_table$genus[name_table$family == family])
    target_node <- lookup_target(anchor_tips[anchor_genera %in% family_genera])
    if (!is.null(target_node)) {
      placement_level  <- "family"
      placement_target <- family
    }
  }
  if (is.null(target_node) && known(order)) {
    order_genera <- unique(name_table$genus[name_table$order == order])
    target_node <- lookup_target(anchor_tips[anchor_genera %in% order_genera])
    if (!is.null(target_node)) {
      placement_level  <- "order"
      placement_target <- order
    }
  }
  if (is.null(target_node)) target_node <- ape::Ntip(tree) + 1L

  tree_before <- tree
  tryCatch({
    tree_height   <- max(phytools::nodeHeights(tree))
    fossil_height <- tree_height - age_ma
    node_height   <- phytools::nodeheight(tree, target_node)
    edge_length   <- fossil_height - node_height

    if (edge_length >= 0) {
      tree <- phytools::bind.tip(
        tree, tip.label = tip_label, where = target_node,
        edge.length = edge_length
      )
    } else {
      # The fossil predates its placement clade: walk up to the branch that
      # was alive at the fossil's age and split it there.
      age_fallback <- "ancestral_branch"
      node  <- target_node
      found <- FALSE
      repeat {
        edge_idx <- which(tree$edge[, 2] == node)
        if (!length(edge_idx)) break
        parent_node   <- tree$edge[edge_idx, 1]
        parent_height <- phytools::nodeheight(tree, parent_node)
        if (parent_height <= fossil_height) {
          edge_limit <- tree$edge.length[edge_idx[1]]
          position   <- phytools::nodeheight(tree, node) - fossil_height
          tolerance  <- sqrt(.Machine$double.eps) * max(1, abs(edge_limit))
          if (position < -tolerance || position > edge_limit + tolerance) {
            stop("computed fossil position lies outside its parent edge")
          }
          position <- max(0, min(position, edge_limit))
          tree <- phytools::bind.tip(
            tree, tip.label = tip_label, where = node, position = position,
            edge.length = 0
          )
          found <- TRUE
          break
        }
        node <- parent_node
      }
      if (!found) {
        stop("fossil age predates the scaffold root; no branch exists at its age")
      }
    }

    placed_height <- phytools::nodeheight(tree, match(tip_label, tree$tip.label))
    tolerance     <- sqrt(.Machine$double.eps) * max(1, abs(fossil_height))
    if (abs(placed_height - fossil_height) > tolerance) {
      stop("grafted tip depth does not match scaffold height minus fossil age")
    }

    list(tree = tree, placed = TRUE, placement_level = placement_level,
         placement_target = placement_target, age_fallback = age_fallback,
         error = NA_character_)
  }, error = function(e) {
    res <- failure(conditionMessage(e))
    res$tree <- tree_before
    res
  })
}
