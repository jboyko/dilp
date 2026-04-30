# Trait columns (dilp-style names) used in the PGLS model.
# Must match the columns computed by dilp_processing().
.pgls_trait_cols <- c(
  "petiole_metric", "ln_leaf_area", "fdr", "margin",
  "perimeter_ratio", "tc_p", "tc_ip", "avg_ta",
  "ta_ba", "ta_p", "ta_ip", "tc_ba"
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
#' fossil species' prediction toward the climate values of its closest extant
#' relatives on the tree.
#'
#' @param specimen_data A data frame containing specimen-level leaf physiognomic
#'   data. Must include all columns required by [dilp_processing()], plus:
#'   * `species` — binomial species name used for phylogenetic placement
#'   * `genus`, `family`, `order` — taxonomic ranks used as placement fallbacks
#'   * `age_ma` — fossil age in millions of years, used for tip placement depth
#'
#'   Placement proceeds genus → family → order → root. A warning is issued for
#'   any species that cannot be placed more precisely than the root.
#'
#' @return A list with five elements:
#' * `results` — a data frame of site-level MAT and MAP estimates with columns
#'   `site`, `MAT.PIP` (°C), `MAP.PIP` (cm), `n_species`, `MAT.PIP.error`
#'   (±°C), `MAP.PIP.error.plus`, and `MAP.PIP.error.minus` (cm). Uncertainty
#'   columns are RMSE-based, from 10-fold LOSO cross-validation in Boyko et al.
#'   (in prep), and do not reflect phylogenetic placement uncertainty.
#' * `species_predictions` — a data frame of per-species-per-site estimates with
#'   columns `site`, `species`, `age_ma`, `MAT.PIP`, and `MAP.PIP`.
#' * `placement_log` — a data frame recording how each fossil species was placed
#'   on the phylogeny (`species`, `placement_level`, `target`).
#' * `processed_leaf_data` — the full processed specimen-level data frame from
#'   [dilp_processing()].
#' * `processed_site_data` — site-level means of leaf physiognomic variables,
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
#' # Add required taxonomy columns to McAbeeExample
#' set.seed(42)
#' fossil_data <- McAbeeExample
#' n <- nrow(fossil_data)
#' fossil_data$species <- paste0("Betula_sp", sample(1:5, n, replace = TRUE))
#' fossil_data$genus   <- "Betula"
#' fossil_data$family  <- "Betulaceae"
#' fossil_data$order   <- "Fagales"
#' fossil_data$age_ma  <- 52.0
#'
#' results <- dilp_pgls(fossil_data)
#' results$results
#' results$species_predictions
#' results$placement_log
#'
#' # Plot sites on a Whittaker biome diagram
#' dilp_whittaker(results)
#'
#' # Check fossil site physiognomy against the calibration space
#' dilp_cca(results)
#' }
dilp_pgls <- function(specimen_data) {

  # ------------------------------------------------------------------
  # 1. Standardise column names
  # ------------------------------------------------------------------
  colnames(specimen_data) <- colnameClean(specimen_data)

  required_tax <- c("species", "genus", "family", "order", "age_ma")
  missing_tax  <- required_tax[!required_tax %in% colnames(specimen_data)]
  if (length(missing_tax) > 0) {
    stop(paste("specimen_data is missing required columns:",
               paste(missing_tax, collapse = ", ")))
  }

  # ------------------------------------------------------------------
  # 2. Process leaf traits
  # ------------------------------------------------------------------
  processed <- dilp_processing(specimen_data)

  # dilp_processing() drops no columns, so taxonomy survives

  # ------------------------------------------------------------------
  # 2b. Fill tooth traits for untoothed specimens
  # For specimens with margin == 1 (entire), tooth-count-derived traits
  # are NA. Fill with 0 (or 1 for perimeter_ratio) so entirely untoothed
  # species aggregate to 0 rather than NaN, matching the training data.
  # ------------------------------------------------------------------
  tooth_cols_zero <- c("tc_ip", "tc_p", "avg_ta", "ta_ba", "ta_p", "ta_ip", "tc_ba")
  untoothed       <- !is.na(processed$margin) & processed$margin == 1

  for (col in tooth_cols_zero) {
    if (col %in% colnames(processed)) {
      processed[[col]] <- ifelse(untoothed & is.na(processed[[col]]),
                                 0, processed[[col]])
    }
  }
  # perimeter_ratio is 1 for untoothed leaves (smooth perimeter; raw = internal)
  if ("perimeter_ratio" %in% colnames(processed)) {
    processed$perimeter_ratio <- ifelse(untoothed & is.na(processed$perimeter_ratio),
                                        1, processed$perimeter_ratio)
  }

  # ------------------------------------------------------------------
  # 3a. Build processed_site_data for compatibility with dilp_cca()
  #     Mirrors dilp() aggregation: specimen -> morphotype -> site
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

  # ------------------------------------------------------------------
  # 3. Aggregate to species grand means (for tree placement)
  #    and species-within-site means (for predictions — PIP+site)
  # ------------------------------------------------------------------
  sp_grand <- processed %>%
    dplyr::group_by(.data$species, .data$genus, .data$family, .data$order) %>%
    dplyr::summarize(
      age_ma = mean(.data$age_ma, na.rm = TRUE),
      dplyr::across(dplyr::all_of(.pgls_trait_cols),
                    \(x) mean(x, na.rm = TRUE)),
      .groups = "drop"
    )

  sp_site <- processed %>%
    dplyr::group_by(.data$site, .data$species) %>%
    dplyr::summarize(
      age_ma = mean(.data$age_ma, na.rm = TRUE),
      dplyr::across(dplyr::all_of(.pgls_trait_cols),
                    \(x) mean(x, na.rm = TRUE)),
      .groups = "drop"
    )

  # ------------------------------------------------------------------
  # 4. Warn about predictors that are still NA after aggregation
  # ------------------------------------------------------------------
  na_counts <- colSums(is.na(sp_site[, .pgls_trait_cols]))
  na_preds  <- na_counts[na_counts > 0]
  if (length(na_preds) > 0) {
    warning(
      "Some species-site observations are missing predictors after aggregation ",
      "and will produce NA predictions.\n  Missing: ",
      paste(paste0(names(na_preds), " (", na_preds, " NA)"), collapse = ", "),
      call. = FALSE
    )
  }

  # ------------------------------------------------------------------
  # 5. Load scaffold tree and place fossil species
  # ------------------------------------------------------------------
  scaffold_path <- system.file("extdata", "tre_scaffold.tre", package = "dilp")
  if (!nzchar(scaffold_path)) {
    stop("Scaffold tree not found. Please reinstall the dilp package.")
  }
  tree <- ape::read.tree(scaffold_path)

  placement   <- .pgls_place_fossils(tree, sp_grand, pgls_name_table)
  tree        <- placement$tree
  placement_log <- placement$log

  fossil_species <- sp_grand$species[sp_grand$species %in% tree$tip.label]
  sp_site        <- sp_site[sp_site$species %in% fossil_species, ]

  if (nrow(sp_site) == 0) {
    stop("No fossil species could be placed on the phylogeny.")
  }

  # ------------------------------------------------------------------
  # 6. Compute phylogenetic adjustment (t(V_cross) %*% Ke)
  # ------------------------------------------------------------------
  idx_extant <- names(pgls_Ke_mat)
  tree_small <- ape::keep.tip(tree, c(idx_extant, fossil_species))
  phylomat   <- ape::vcv(tree_small)

  V_cross <- phylomat[idx_extant, fossil_species, drop = FALSE]

  phylo_adj_mat <- stats::setNames(
    as.numeric(t(V_cross) %*% pgls_Ke_mat[idx_extant]),
    fossil_species
  )
  phylo_adj_map <- stats::setNames(
    as.numeric(t(V_cross) %*% pgls_Ke_map[idx_extant]),
    fossil_species
  )

  # ------------------------------------------------------------------
  # 7. Build design matrix and predict per species-site
  # ------------------------------------------------------------------
  X <- cbind(
    `(Intercept)` = 1,
    as.matrix(sp_site[, .pgls_trait_cols])
  )
  X <- X[, names(pgls_beta_mat), drop = FALSE]

  yhat_mat <- as.numeric(X %*% pgls_beta_mat) +
                phylo_adj_mat[sp_site$species]
  yhat_map <- as.numeric(X %*% pgls_beta_map) +
                phylo_adj_map[sp_site$species]

  sp_site$MAT.PIP <- yhat_mat
  sp_site$MAP.PIP <- exp(yhat_map)

  # ------------------------------------------------------------------
  # 8. Aggregate species predictions to site level
  # ------------------------------------------------------------------
  results <- sp_site %>%
    dplyr::group_by(.data$site) %>%
    dplyr::summarize(
      MAT.PIP   = mean(.data$MAT.PIP,  na.rm = TRUE),
      MAP.PIP   = mean(.data$MAP.PIP,  na.rm = TRUE),
      n_species = dplyr::n(),
      .groups   = "drop"
    )

  # ------------------------------------------------------------------
  # 9. Attach RMSE-based uncertainty from 10-fold LOSO cross-validation
  #    (Boyko et al., in prep). MAT error is symmetric (°C); MAP errors
  #    are asymmetric, back-transformed from the log scale.
  # ------------------------------------------------------------------
  log_map_pip <- log(results$MAP.PIP)
  results$MAT.PIP.error       <- 3.42
  results$MAP.PIP.error.plus  <- exp(log_map_pip + 0.52) - results$MAP.PIP
  results$MAP.PIP.error.minus <- results$MAP.PIP - exp(log_map_pip - 0.52)

  message(
    "Uncertainty columns (MAT.PIP.error, MAP.PIP.error.*) are RMSE-based, ",
    "derived from 10-fold leave-one-site-out cross-validation reported in ",
    "Boyko et al. (in prep). They do not account for phylogenetic placement uncertainty."
  )

  list(
    results             = as.data.frame(results),
    species_predictions = as.data.frame(
      sp_site[, c("site", "species", "age_ma", "MAT.PIP", "MAP.PIP")]
    ),
    placement_log       = placement_log,
    processed_leaf_data = processed,
    processed_site_data = as.data.frame(processed_site_data)
  )
}


# ----------------------------------------------------------------------
# Internal: graft fossil species onto scaffold tree
# ----------------------------------------------------------------------
.pgls_place_fossils <- function(tree, sp_grand, name_table) {

  h          <- max(phytools::nodeHeights(tree))
  tip_genera <- sapply(strsplit(tree$tip.label, "_"), `[`, 1)
  log_rows   <- vector("list", nrow(sp_grand))

  for (j in seq_len(nrow(sp_grand))) {
    fossil_name <- sp_grand$species[j]
    age_ma      <- sp_grand$age_ma[j]

    if (fossil_name %in% tree$tip.label) {
      log_rows[[j]] <- data.frame(species       = fossil_name,
                                  placement_level = "already_present",
                                  target         = fossil_name,
                                  stringsAsFactors = FALSE)
      next
    }

    target_node    <- NULL
    placement_level <- NULL
    placement_target <- NULL

    # Genus match
    genus_tips <- tree$tip.label[tip_genera == sp_grand$genus[j]]
    if (length(genus_tips) >= 2) {
      target_node     <- ape::getMRCA(tree, genus_tips)
      placement_level <- "genus"
      placement_target <- sp_grand$genus[j]
    } else if (length(genus_tips) == 1) {
      target_node     <- which(tree$tip.label == genus_tips[1])
      placement_level <- "genus"
      placement_target <- sp_grand$genus[j]
    }

    # Family fallback
    if (is.null(target_node) &&
        nzchar(sp_grand$family[j]) && sp_grand$family[j] != "unknown") {
      fam_genera <- unique(name_table$genus[name_table$family == sp_grand$family[j]])
      fam_tips   <- tree$tip.label[tip_genera %in% fam_genera]
      if (length(fam_tips) >= 2) {
        target_node     <- ape::getMRCA(tree, fam_tips)
        placement_level <- "family"
        placement_target <- sp_grand$family[j]
      } else if (length(fam_tips) == 1) {
        target_node     <- which(tree$tip.label == fam_tips[1])
        placement_level <- "family"
        placement_target <- sp_grand$family[j]
      }
    }

    # Order fallback
    if (is.null(target_node) &&
        nzchar(sp_grand$order[j]) && sp_grand$order[j] != "unknown") {
      ord_genera <- unique(name_table$genus[name_table$order == sp_grand$order[j]])
      ord_tips   <- tree$tip.label[tip_genera %in% ord_genera]
      if (length(ord_tips) >= 2) {
        target_node     <- ape::getMRCA(tree, ord_tips)
        placement_level <- "order"
        placement_target <- sp_grand$order[j]
      } else if (length(ord_tips) == 1) {
        target_node     <- which(tree$tip.label == ord_tips[1])
        placement_level <- "order"
        placement_target <- sp_grand$order[j]
      }
    }

    # Root fallback
    if (is.null(target_node)) {
      warning("No taxonomy match for '", fossil_name,
              "'. Placing at root.", call. = FALSE)
      target_node     <- length(tree$tip.label) + 1L
      placement_level <- "root"
      placement_target <- "root"
    }

    log_rows[[j]] <- data.frame(species        = fossil_name,
                                placement_level = placement_level,
                                target         = placement_target,
                                stringsAsFactors = FALSE)

    node_ht  <- phytools::nodeheight(tree, target_node)
    edge_len <- (h - age_ma) - node_ht

    if (edge_len >= 0) {
      tree <- phytools::bind.tip(tree, tip.label = fossil_name,
                                 where = target_node,
                                 edge.length = edge_len)
    } else {
      # Fossil age predates target node — walk up to find a spanning branch
      fossil_ht <- h - age_ma
      node      <- target_node
      placed    <- FALSE

      repeat {
        parent_idx  <- which(tree$edge[, 2] == node)
        if (length(parent_idx) == 0) break
        parent_node <- tree$edge[parent_idx, 1]
        parent_ht   <- phytools::nodeheight(tree, parent_node)
        if (parent_ht <= fossil_ht) {
          position <- phytools::nodeheight(tree, node) - fossil_ht
          tree <- phytools::bind.tip(tree, tip.label = fossil_name,
                                     where = node, position = position,
                                     edge.length = 0)
          placed <- TRUE
          break
        }
        node <- parent_node
      }

      if (!placed) {
        warning("No spanning branch for '", fossil_name,
                "'. Placing at root.", call. = FALSE)
        tree <- phytools::bind.tip(tree, tip.label = fossil_name,
                                   where = length(tree$tip.label) + 1L,
                                   edge.length = 0.001)
      }
    }

    tip_genera <- sapply(strsplit(tree$tip.label, "_"), `[`, 1)
  }

  log_df <- do.call(rbind, log_rows[!sapply(log_rows, is.null)])
  list(tree = tree, log = log_df)
}
