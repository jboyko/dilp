#' Process raw leaf physiognomic data
#'
#' @description
#' `dilp_processing()` will typically only be called internally by `dilp()`.
#' However, it can be used on its own to generate and view a processed DiLP
#' dataset that includes raw and derived physiognomic values useful for DiLP and
#' other physiognomic analyses.  Returns a data frame.
#'
#' @param specimen_data A data frame containing specimen level leaf physiognomic
#' data. A good reference for how to put together the data: \code{\link{McAbeeExample}}
#' @return A data frame containing cleaned and processed specimen level leaf
#' physiognomic data.  New variables calculated are:
#' * Leaf area
#' * Feret diameter
#' * Feret diameter ratio (FDR)
#' * Raw blade perimeter corrected (Raw blade perimeter - length of cut perimeter)
#' * Internal raw blade perimeter corrected (Internal raw blade perimeter - length of cut perimeter)
#' * Total tooth count
#' * Total tooth count : internal perimeter (TC:IP)
#' * Perimeter ratio
#' * Petiole metric
#' * Aspect ratio
#' * Shape factor
#' * Compactness
#' * Tooth area
#' * Tooth area : perimeter (TA:P)
#' * Tooth area: internal perimeter (TA:IP)
#' * Tooth area : blade area (TA:BA)
#' * Average primary tooth area (Avg TA)
#' * Tooth count : blade area (TC:BA)
#' * Tooth count : perimeter (TC:P)
#' @export
#'
#' @examples
#' dilp_dataset <- dilp_processing(McAbeeExample)
#' dilp_dataset
dilp_processing <- function(specimen_data) {
  specimen_data <- dplyr::as_tibble(specimen_data)
  colnames(specimen_data) <- colnameClean(specimen_data)

  required_columns <- c(
    "site", "specimen_number", "morphotype", "margin", "feret", "blade_area",
    "raw_blade_perimeter", "internal_raw_blade_perimeter", "length_of_cut_perimeter",
    "no_of_primary_teeth", "no_of_subsidiary_teeth"
  )

  recommended_columns <- c(
    "petiole_width", "petiole_area", "blade_perimeter",
    "minimum_feret", "raw_blade_area", "internal_raw_blade_area"
  )
  missing_columns <- required_columns[!required_columns %in% colnames(specimen_data)]
  if (length(missing_columns) > 0) {
    stop(paste("specimen_data is missing required columns:", stringr::str_flatten(missing_columns, collapse = ", ")))
  } else {
    other_columns <- recommended_columns[!recommended_columns %in% colnames(specimen_data)]
    if (length(other_columns) > 0) {
      warning(paste("specimen_data is missing recommended columns and has filled them with NAs:", stringr::str_flatten(other_columns, collapse = ", ")), call. = FALSE)
      for (column in other_columns) {
        specimen_data[[column]] <- NA
      }
    }

    # Identify 0.5 margin species
    mixed_margins <- specimen_data %>%
      dplyr::group_by(.data$site, .data$morphotype) %>%
      dplyr::filter(1 %in% .data$margin & 0 %in% .data$margin) %>%
      dplyr::filter(.data$margin == 1) %>%
      dplyr::select("site", "morphotype", "specimen_number")


    if (length(mixed_margins$specimen_number) > 0) {
      for (i in 1:length(mixed_margins$specimen_number)) {
        specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$length_of_cut_perimeter <- 0
        specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$no_of_primary_teeth <- 0
        specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$no_of_subsidiary_teeth <- 0
      }
    }


    ## For specimens where blade area was measured, but there is an NA for petiole area, convert NA to 0 to permit the subsequent addition
    temp3 <- specimen_data$petiole_area
    specimen_data$petiole_area <- ifelse(specimen_data$blade_area > 0 & is.na(specimen_data$petiole_area), 0, specimen_data$petiole_area)
    # Sum petiole area and blade area
    specimen_data$leaf_area <- specimen_data$blade_area + specimen_data$petiole_area # if petiole area was measured but there is a NA for blade area, this addition will appropriately return a NA for leaf area

    # Feret Diameter Ratio
    specimen_data$feret_diameter <- 2 * sqrt(specimen_data$leaf_area / pi)
    specimen_data$fdr <- specimen_data$feret_diameter / specimen_data$feret

    # change NA's to zeros for length.of.cut.perimeter so that derived variables can be calculated
    temp1 <- specimen_data$length_of_cut_perimeter
    specimen_data$length_of_cut_perimeter <- ifelse(is.na(specimen_data$length_of_cut_perimeter), 0, specimen_data$length_of_cut_perimeter)

    # Corrected perimeters
    specimen_data$raw_blade_perimeter <- ifelse(is.na(specimen_data$raw_blade_perimeter) & !is.na(specimen_data$blade_perimeter) & specimen_data$margin == 0 & specimen_data$length_of_cut_perimeter == 0,
      specimen_data$blade_perimeter, specimen_data$raw_blade_perimeter
    )
    specimen_data$raw_blade_perimeter_corrected <- specimen_data$raw_blade_perimeter - specimen_data$length_of_cut_perimeter
    specimen_data$internal_raw_blade_perimeter_corrected <- specimen_data$internal_raw_blade_perimeter - specimen_data$length_of_cut_perimeter

    # Corrected raw blade area if there is no cut perimeter on a toothed leaf.
    specimen_data$raw_blade_area <- ifelse(is.na(specimen_data$raw_blade_area) & specimen_data$margin == 0 & specimen_data$length_of_cut_perimeter == 0,
      specimen_data$blade_area, specimen_data$raw_blade_area
    )

    ## Toothed variables
    # Total tooth count
    temp2 <- specimen_data$no_of_subsidiary_teeth
    specimen_data$no_of_subsidiary_teeth[is.na(specimen_data$no_of_subsidiary_teeth)] <- 0 # if it is left as NA then total tooth count will also return NA
    specimen_data$total_tooth_count <- ifelse(specimen_data$margin == 1, NA, specimen_data$no_of_primary_teeth + specimen_data$no_of_subsidiary_teeth)
    if (length(mixed_margins$specimen_number) > 0) {
      for (i in 1:length(mixed_margins$specimen_number)) {
        specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$total_tooth_count <- 0
      }
    }

    # Total tooth count : internal perimeter
    specimen_data$tc_ip <- specimen_data$total_tooth_count / specimen_data$internal_raw_blade_perimeter_corrected
    if (length(mixed_margins$specimen_number) > 0) {
      for (i in 1:length(mixed_margins$specimen_number)) {
        specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$tc_ip <- 0
      }
    }

    # Perimeter ratio
    specimen_data$perimeter_ratio <- specimen_data$raw_blade_perimeter_corrected / specimen_data$internal_raw_blade_perimeter_corrected
    if (length(mixed_margins$specimen_number) > 0) {
      for (i in 1:length(mixed_margins$specimen_number)) {
        specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$perimeter_ratio <- 1
      }
    }

    #### Apply required natural logs to appropriate variables
    specimen_data$ln_leaf_area <- log(100 * specimen_data$leaf_area) # Leaf area is expressed in mm2 in DiLP MLR and SLR models, so here also converting from cm2 to mm2
    specimen_data$ln_pr <- log(specimen_data$perimeter_ratio)
    specimen_data$ln_tc_ip <- log(specimen_data$tc_ip)
    if (length(mixed_margins$specimen_number) > 0) {
      for (i in 1:length(mixed_margins$specimen_number)) {
        specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$ln_tc_ip <- NA
      }
    }

    # Petiole metric for reconstructing leaf mass per area
    specimen_data$petiole_metric <- (specimen_data$petiole_width^2) / specimen_data$leaf_area
    # Aspect ratio
    specimen_data$aspect_ratio <- specimen_data$feret / specimen_data$minimum_feret
    # Shape factor
    specimen_data$shape_factor <- 4 * pi * (specimen_data$blade_area / (specimen_data$blade_perimeter^2))
    # Compactness
    specimen_data$compactness <- (specimen_data$blade_perimeter^2) / specimen_data$blade_area
    # Tooth area
    specimen_data$tooth_area <- specimen_data$raw_blade_area - specimen_data$internal_raw_blade_area
    if (length(mixed_margins$specimen_number) > 0) {
      for (i in 1:length(mixed_margins$specimen_number)) {
        specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$tooth_area <- 0
      }
    }

    # Tooth area : perimeter
    specimen_data$ta_p <- specimen_data$tooth_area / specimen_data$raw_blade_perimeter_corrected
    if (length(mixed_margins$specimen_number) > 0) {
      for (i in 1:length(mixed_margins$specimen_number)) {
        specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$ta_p <- 0
      }
    }

    # Tooth area : internal perimeter
    specimen_data$ta_ip <- specimen_data$tooth_area / specimen_data$internal_raw_blade_perimeter_corrected
    if (length(mixed_margins$specimen_number) > 0) {
      for (i in 1:length(mixed_margins$specimen_number)) {
        specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$ta_ip <- 0
      }
    }

    # Tooth area : blade area
    specimen_data$ta_ba <- specimen_data$tooth_area / specimen_data$raw_blade_area
    if (length(mixed_margins$specimen_number) > 0) {
      for (i in 1:length(mixed_margins$specimen_number)) {
        specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$ta_ba <- 0
      }
    }

    # Average primary tooth area
    specimen_data$avg_ta <- specimen_data$tooth_area / specimen_data$no_of_primary_teeth # double check
    if (length(mixed_margins$specimen_number) > 0) {
      for (i in 1:length(mixed_margins$specimen_number)) {
        specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$avg_ta <- 0
      }
    }

    # Tooth count : blade area
    specimen_data$tc_ba <- specimen_data$total_tooth_count / specimen_data$raw_blade_area
    if (length(mixed_margins$specimen_number) > 0) {
      for (i in 1:length(mixed_margins$specimen_number)) {
        specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$tc_ba <- 0
      }
    }

    # tooth count : perimeter
    specimen_data$tc_p <- specimen_data$total_tooth_count / specimen_data$raw_blade_perimeter_corrected
    if (length(mixed_margins$specimen_number) > 0) {
      for (i in 1:length(mixed_margins$specimen_number)) {
        specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$tc_p <- 0
      }
    }


    # Return length of cut perimeter and no subsidiary teeth back to original values to keep averages intact
    specimen_data$length_of_cut_perimeter <- temp1
    specimen_data$no_of_subsidiary_teeth <- temp2
    specimen_data$petiole_area <- temp3




    return(specimen_data)
  }
}

#' Check for common errors in DiLP measurements
#'
#' @description
#' `dilp_errors()` will typically only be called internally by [dilp()].
#' However, it can be used on its own to evaluate errors that commonly occur
#' during the data collection and processing steps.  A `dilp_errors()` call
#' will nearly always follow a [dilp_processing()] call.  Returns a data frame.
#'
#' @param specimen_data Processed specimen level leaf physiognomic data.  The
#' structure should match the structure of the output from [dilp_processing()]
#'
#' @return A 7 by X data frame.  Each row shows a common error, and which specimens
#' from the input dataset are tripping it.
#' @export
#'
#' @examples
#' # Check for errors in the provided McAbeeExample dataset.
#' dilp_dataset <- dilp_processing(McAbeeExample)
#' dilp_errors <- dilp_errors(dilp_dataset)
#' dilp_errors
dilp_errors <- function(specimen_data) {
  ### This data frame will be filled with the results of the following error checks
  errors <- data.frame()
  #Flag specimens in mixed margin morphotypes
  mixed_margins <- specimen_data %>%
    dplyr::group_by(.data$site, .data$morphotype) %>%
    dplyr::filter(1 %in% .data$margin & 0 %in% .data$margin) %>%
    dplyr::filter(.data$margin == 1) %>%
    dplyr::select("site", "morphotype", "specimen_number")

  #Isolate entire margin specimens
  check1 <- specimen_data %>% dplyr::filter(.data$margin == 1)

  #Originally we were checking total tooth count, but it is automatically converted to NA in the processing step in case people put 0's in. Thus here we instead check primary and subsidiary teeth separately. Presumably this is needed for mixed margins.
  error <- specimen_data[which(check1$no_of_primary_teeth > -1 & !(check1$specimen_number %in% mixed_margins$specimen_number)), "specimen_number"] #The mixed margin bit ensures that specimens in a mixed margin morphotype retain their tooth count and area values of 0, and no error is flagged for it
  temp1 <- data.frame(Check = "Entire tooth count not NA", specimen_number = ifelse(nrow(error) != 0, paste0(error$specimen_number, collapse = ", "), "No errors found"))
  error <- specimen_data[which(check1$no_of_subsidiary_teeth > -1 & !(check1$specimen_number %in% mixed_margins$specimen_number)), "specimen_number"] #The mixed margin bit ensures that specimens in a mixed margin morphotype retain their tooth count and area values of 0, and no error is flagged for it
  temp1.1 <- data.frame(Check = "Entire tooth count not NA", specimen_number = ifelse(nrow(error) != 0, paste0(error$specimen_number, collapse = ", "), "No errors found"))

  colnames(temp1) <- c("Check", "specimen_number")
  colnames(temp1.1) <- c("Check", "specimen_number")
  temp1 <- rbind(temp1, temp1.1)
  #If the number of unique values in specimen number is greater than 1, errors were reported
  #If it was only reported for primary or subsidiary, one of the returns will still say "No errors found", which needs to be deleted
  if(length(unique(temp1$specimen_number))>1){temp1 <- dplyr::filter(temp1, .data$specimen_number != "No errors found")}
  #At the end, keep only unique specimen numbers, to ensure either specimen numbers of "No errors found" is not repeated
  temp1 <- dplyr::distinct(temp1)

  error <- specimen_data[which(check1$tc_ip > -1 & !(check1$specimen_number %in% mixed_margins$specimen_number)), "specimen_number"]
  temp2 <- data.frame(Check = "Entire tooth count : IP not NA", specimen_number = ifelse(nrow(error) != 0, paste0(error$specimen_number, collapse = ", "), "No errors found"))
  colnames(temp2) <- c("Check", "specimen_number")

  error <- data.frame("specimen_number" = check1$specimen_number[which(check1$perimeter_ratio > -1 & !(check1$specimen_number %in% mixed_margins$specimen_number))])

  temp3 <- data.frame(Check = "Entire perimeter ratio not NA", specimen_number = ifelse(nrow(error) != 0, paste0(error$specimen_number, collapse = ", "), "No errors found"))
  colnames(temp3) <- c("Check", "specimen_number")

  check2 <- tidyr::drop_na(specimen_data, "fdr")

  error <- specimen_data[which(check2$fdr < 0 | check2$fdr > 1), "specimen_number"]
  temp4 <- data.frame(Check = "FDR not between 0-1", specimen_number = ifelse(nrow(error) != 0, paste0(error$specimen_number, collapse = ", "), "No errors found"))
  colnames(temp4) <- c("Check", "specimen_number")

  error <- specimen_data[which(specimen_data$internal_raw_blade_perimeter_corrected > specimen_data$raw_blade_perimeter_corrected), "specimen_number"]
  temp5 <- data.frame(Check = "External perimeter not larger than internal perimeter", specimen_number = ifelse(nrow(error) != 0, paste0(error$specimen_number, collapse = ", "), "No errors found"))
  colnames(temp5) <- c("Check", "specimen_number")

  error <- specimen_data[which(specimen_data$minimum_feret > specimen_data$feret), "specimen_number"]
  temp6 <- data.frame(Check = "Feret is not larger than minimum Feret", specimen_number = ifelse(nrow(error) != 0, paste0(error$specimen_number, collapse = ", "), "No errors found"))
  colnames(temp6) <- c("Check", "specimen_number")

  error <- as.data.frame(specimen_data[which(specimen_data$perimeter_ratio < 1 & !(specimen_data$specimen_number %in% mixed_margins$specimen_number)), "specimen_number"])
  temp7 <- data.frame(Check = "Perimeter ratio less than 1", specimen_number = ifelse(nrow(error) != 0, paste0(error$specimen_number, collapse = ", "), "No errors found"))
  colnames(temp7) <- c("Check", "specimen_number")

  errors <- dplyr::bind_rows(temp1, temp2, temp3, temp4, temp5, temp6, temp7)
  rownames(errors) <- NULL

  if (length(unique(errors$specimen_number)) > 1){
    warning("Errors found in dataset. Check $errors for list.")
  }
  return(errors)
}

#' Identify outlier specimens
#'
#' @description
#' `dilp_outliers()` is called internally by [dilp()].
#' However, it can be used on its own to flag specimens that may have been
#' reported, measured, or prepared incorrectly.  `dilp_outliers()` returns a data frame
#' listing specimens that have unusually high or low values for the four key
#' parameters used in DiLP analyses. This includes whether a specimen is an outlier
#' in the entire dataset, or among other specimens in the same morphotype. If flagged, we suggest looking at the
#' raw measurements and prepped specimen and evaluating if the data is in error or is correct. If in error, the specimen will
#' need to be reprepared and/or remeasured, and the updated datasheet re-read back into R.
#'
#'
#' @param specimen_data Processed specimen level leaf physiognomic data.  The
#' structure should match the structure of the output from [dilp_processing()]
#'
#' @return A 4 by X data frame. Each row represents one of the DiLP parameters,
#' and the specimens that are outliers for that parameter.
#' @export
#'
#' @examples
#' # Check for outliers in the provided McAbeeExample dataset. Each
#' # of these outliers has been manually re-examined and was found acceptable.
#' dilp_dataset <- dilp_processing(McAbeeExample)
#' dilp_outliers <- dilp_outliers(dilp_dataset)
#' dilp_outliers
dilp_outliers <- function(specimen_data) {
  vars <- c("fdr", "tc_ip", "leaf_area", "perimeter_ratio") # DiLP variables
  outliers <- data.frame()

  #####Outliers within entire dataset
  for (i in 1:length(vars)) {
    temp <- specimen_data
    colnames(temp)[colnames(temp) == vars[i]] <- "trait" # rename variable of focus vars[i]
    temp.outliers <- grDevices::boxplot.stats(temp$trait)$out # check it for outliers
    temp.specimen <- temp[which(temp$trait %in% c(temp.outliers)), c("site", "specimen_number", "morphotype")] # determine specimen numbers for any outliers
    temp.output <- as.data.frame(temp.specimen)
    if(nrow(temp.output)>0){temp.output$outlier.entire.dataset <- vars[i]}
    outliers <- if(nrow(temp.output)>0){dplyr::bind_rows(outliers, temp.output)} # bind to summary table
  }

  #####Outliers within morphotype
  morphs <- unique(specimen_data$morphotype)
  outliers.2 <- data.frame()
  for (j in 1:length(morphs)) {
    temp <- specimen_data
    temp.morph <- dplyr::filter(temp, .data$morphotype == morphs[j])

    for (i in 1:length(vars)) {
      temp.morph2 <- temp.morph
      colnames(temp.morph2)[colnames(temp.morph2) == vars[i]] <- "trait" # rename variable of focus vars[i]
      temp.outliers <- grDevices::boxplot.stats(temp.morph2$trait)$out # check it for outliers
      temp.specimen <- temp.morph2[which(temp.morph2$trait %in% c(temp.outliers)), c("site", "specimen_number", "morphotype")] # determine specimen numbers for any outliers
      temp.output <- as.data.frame(temp.specimen)
      if(nrow(temp.output)>0){temp.output$outlier.morphotype <- vars[i]}
      if(nrow(temp.output)>0){outliers.2 <- dplyr::bind_rows(outliers.2, temp.output)} # bind to summary table
    }
  }
  #Merge entire dataset and morphotye datasets
  outliers <- merge(outliers, outliers.2, by=c("site", "specimen_number", "morphotype"), all=TRUE)
  #Replace NA's with "No outliers"
  outliers$outlier.entire.dataset[is.na(outliers$outlier.entire.dataset)] <- "No outliers"
  outliers$outlier.morphotype[is.na(outliers$outlier.morphotype)] <- "No outliers"
  ##Sort by specimen number
  if(length(outliers) != 0) {
    warning("Outliers found. Please evaluate $outliers for possible wrong measurements")
  }
  return(outliers)

}

#' Generate DiLP results
#'
#' @description
#' `dilp()` processes raw leaf physiognomic data, checks for common
#' errors/outliers, and returns the processed data, keys to finding potential
#' errors or outliers, and paleoclimate reconstructions.
#'
#' @param specimen_data A data frame containing specimen level leaf physiognomic
#' data. See Lowe et al. 2024 for more information on how to collect this data.
#' A good reference for how to put together the data: \code{\link{McAbeeExample}}
#'
#' Required columns:
#'    * site
#'    * specimen_number
#'    * morphotype
#'    * margin
#'    * feret
#'    * blade_area
#'    * raw_blade_perimeter
#'    * internal_raw_blade_perimeter
#'    * length_of_cut_perimeter
#'    * no_of_primary_teeth
#'    * no_of_subsidiary_teeth
#'
#'  Recommended columns:
#'    * petiole_width
#'    * petiole_area
#'    * blade_perimeter
#'    * minimum_feret
#'    * raw_blade_area
#'    * internal_raw_blade_area
#'
#' @param params Either a string referring to one of two preloaded parameter sets
#' of a list of custom parameters (same format as the list below).
#'
#' Preloaded parameter sets are "PeppeGlobal" and "PeppeNH" which are calibrated based on
#' global and northern hemisphere data respectively. Allen et al. (2020) illustrates a situation
#' in which the northern hemisphere parameters may be preferable.  The "PeppeNH" parameters
#' only estimate MAT.  Use "PeppeGlobal" for all MAP estimates. Defaults to "PeppeGlobal" as follows (Peppe et al. 2011):
#'
#'    * MAT.MLR.M = 0.21,
#'    * MAT.MLR.FDR = 42.296,
#'    * MAT.MLR.TC.IP = -2.609,
#'    * MAT.MLR.constant = -16.004,
#'    * MAT.MLR.error = 4,
#'    * MAT.SLR.M = 0.204,
#'    * MAT.SLR.constant = 4.6,
#'    * MAT.SLR.error = 4.9,
#'    * MAP.MLR.LA = 0.298,
#'    * MAP.MLR.TC.IP = 0.279,
#'    * MAP.MLR.PR = -2.717,
#'    * MAP.MLR.constant = 3.033,
#'    * MAP.MLR.SE = 0.6,
#'    * MAP.SLR.LA = 0.283,
#'    * MAP.SLR.constant = 2.92,
#'    * MAP.SLR.SE = 0.61
#'
#' @param subsite_cols A vector or list of columns present in `specimen_data` to calculate
#' paleoclimate estimates for.  A completely optional parameter - allows different groupings of
#' specimens to be tested, or comparisons of paleoclimate estimates at different levels of grouping.
#' Adds additional estimates to $results.
#'
#' @return A list of tables that includes all pertinent DiLP
#' information:
#'
#' * processed_leaf_data: the full set of cleaned and newly calculated leaf
#' physiognomic data that is necessary for DiLP analysis. See [dilp_processing()]
#' for more information.
#' * processed_morphotype_data: morphospecies-site pair means for all leaf
#' physiognomic data.
#' * processed_site_data: site means for all leaf physiognomic data.
#' * errors: lists any specimens that may be causing common errors in DiLP
#' calculations. See [dilp_errors()] for more information.
#' * outliers: flags outliers in variables used for DiLP analysis that may
#' represent incorrect data.  See [dilp_outliers()] for more information.
#' * results: climate reconstructions of MAT and MAP using single and multi-linear
#' regressions.
#'
#' @references
#' * Allen, S. E., Lowe, A. J., Peppe, D. J., & Meyer, H. W. (2020). Paleoclimate and paleoecology of the latest Eocene Florissant flora of central Colorado, USA. Palaeogeography, Palaeoclimatology, Palaeoecology, 551, 109678.
#' * Peppe, D.J., Royer, D.L., Cariglino, B., Oliver, S.Y., Newman, S., Leight, E., Enikolopov, G., Fernandez-Burgos, M., Herrera, F., Adams, J.M., Correa, E., Currano, E.D., Erickson, J.M., Hinojosa, L.F., Hoganson, J.W., Iglesias, A., Jaramillo, C.A., Johnson, K.R., Jordan, G.J., Kraft, N.J.B., Lovelock, E.C., Lusk, C.H., Niinemets, Ü., Peñuelas, J., Rapson, G., Wing, S.L. and Wright, I.J. (2011), Sensitivity of leaf size and shape to climate: global patterns and paleoclimatic applications. New Phytologist, 190: 724-739. https://doi.org/10.1111/j.1469-8137.2010.03615.x
#' * Lowe. A.J., Flynn, A.G., Butrim, M.J., Baumgartner, A., Peppe, D.J., and Royer, D.L. (2024), Reconstructing terrestrial paleoclimate and paleoecology with fossil leaves using Digital Leaf Physiognomy and leaf mass per area.  J. Vis. Exp. (212), e66838, doi:10.3791/66838 (2024).
#' @export
#'
#' @examples
#' dilp_results <- dilp(McAbeeExample)
#' dilp_results$processed_leaf_data
#' dilp_results$processed_morphotype_data
#' dilp_results$processed_site_data
#' dilp_results$errors
#' dilp_results$outliers
#' dilp_results$results
dilp <- function(specimen_data, params = "PeppeGlobal", subsite_cols = NULL) {
  subsite_cols <- subsite_cols %>%
    stringr::str_trim() %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("[.]", " ") %>%
    stringr::str_replace_all("[ ]", "_") %>%
    stringr::str_replace_all("__", "_")
  colnames(specimen_data) <- colnameClean(specimen_data)
  processed_specimen_data <- dilp_processing(specimen_data)
  errors <- dilp_errors(processed_specimen_data)
  outliers <- dilp_outliers(processed_specimen_data)
  if (is.list(params)) {

  } else {
    params <- grab_regression(params, "dilp")
  }

  ####### Morphotype average by site
  dilp_morphotype <- processed_specimen_data %>%
    dplyr::select(-"specimen_number") %>%
    dplyr::summarize(.by = c(.data$site, .data$morphotype), dplyr::across(dplyr::where(~ !is.character(.)), \(x) mean(x, na.rm = TRUE)))

  ##### Morphotypes that have variable leaf margin states require a margin state of 0.5

  # is just 1 and 0 listed? If so, the following finds margin state values between 0 and 1 and replaces them with 0.5
  if (length(unique(stats::na.omit(dilp_morphotype$margin))) > 2) {
    dilp_morphotype$margin[dilp_morphotype$margin > 0 & dilp_morphotype$margin < 1] <- 0.5
    if (length(unique(dilp_morphotype$margin)) > 2) {
      warning("Margin states outside the bounds of [0 - 1] present. Non 0 and 1 values converted to 0.5")
    }
  }

  ####### Site average
  dilp_site <- dilp_morphotype %>%
    dplyr::group_by(.data$site) %>%
    dplyr::select(-"morphotype") %>%
    dplyr::summarize(dplyr::across(dplyr::where(~ !is.character(.)), \(x) mean(x, na.rm = TRUE)))

  ### Convert site margin from proportion to percentage
  dilp_site$margin <- dilp_site$margin * 100

  ### Conditionally set site tc_ip to 0 and perimeter_ratio to 1 if all morphotypes are untoothed
  dilp_site$tc_ip <- ifelse(dilp_site$margin == 100, 0, dilp_site$tc_ip)
  dilp_site$perimeter_ratio <- ifelse(dilp_site$margin == 100, 1, dilp_site$perimeter_ratio)
  # dilp_site$ln_tc_ip <- ifelse(dilp_site$margin == 100, -20.7, dilp_site$ln_tc_ip)
  # dilp_site$ln_pr <- ifelse(dilp_site$margin == 100, 0, dilp_site$ln_pr)

  ## A loop is constructed to handle spreadsheets that contain either multiple sites or a single site. For the latter, the loop will run only once.
  sites <- c(unique(dilp_site$site))
  Results <- data.frame()

  for (i in 1:length(sites)) {
    temp <- dilp_site[dilp_site$site == sites[i], ] # isolate site

    #### MAT
    # MLR
    MAT.MLR <- (temp$margin * params$MAT.MLR.M) + (temp$fdr * params$MAT.MLR.FDR) + (temp$tc_ip * params$MAT.MLR.TC.IP) + params$MAT.MLR.constant

    # SLR
    MAT.SLR <- (temp$margin * params$MAT.SLR.M) + params$MAT.SLR.constant

    #### MAP
    # MLR
    MAP.MLR.exp <- (temp$ln_leaf_area * params$MAP.MLR.LA) + (temp$ln_tc_ip * params$MAP.MLR.TC.IP) + (temp$ln_pr * params$MAP.MLR.PR) + params$MAP.MLR.constant
    MAP.MLR <- exp(MAP.MLR.exp)

    MAP.MLR.error.plus <- (exp(MAP.MLR.exp + params$MAP.MLR.SE)) - MAP.MLR
    MAP.MLR.error.minus <- MAP.MLR - (exp(MAP.MLR.exp - params$MAP.MLR.SE))

    # SLR

    MAP.SLR.exp <- (temp$ln_leaf_area * params$MAP.SLR.LA) + params$MAP.SLR.constant
    MAP.SLR <- exp(MAP.SLR.exp)
    MAP.SLR.error.plus <- (exp(MAP.SLR.exp + params$MAP.SLR.SE)) - MAP.SLR
    MAP.SLR.error.minus <- MAP.SLR - (exp(MAP.SLR.exp - params$MAP.SLR.SE))

    # Results table

    temp_output <- data.frame(
      site = temp$site, margin = temp$margin, fdr = temp$fdr, tc_ip = temp$tc_ip, ln_leaf_area = temp$ln_leaf_area, ln_tc_ip = temp$ln_tc_ip, ln_pr = temp$ln_pr, MAT.MLR, MAT.MLR.error = params$MAT.MLR.error, MAT.SLR, MAT.SLR.error = params$MAT.SLR.error, MAP.MLR,
      MAP.MLR.error.plus, MAP.MLR.error.minus, MAP.SLR, MAP.SLR.error.plus, MAP.SLR.error.minus
    )
    Results <- rbind(Results, temp_output)
  }

  for (subsite in subsite_cols) {
    if (is.null(specimen_data[[subsite]])) {
      stop(paste("Subsite column:", subsite, "not found in specimen_data."))
    } else {
      subsite_data <- specimen_data %>%
        dplyr::select(-.data$site) %>%
        dplyr::mutate(site = .data[[subsite]])
      temp <- dilp(subsite_data)
      Results <- rbind(Results, temp$results)
    }
  }

  return(list(processed_leaf_data = processed_specimen_data, processed_morphotype_data = dilp_morphotype, processed_site_data = dilp_site, errors = errors, outliers = outliers, results = Results))
}
