#' Test if site leaf physiognomy falls within the physiognomic space of the DiLP calibration dataset
#'
#' `dilp_cca` plots a canonical correspondence analysis (CCA) ordination of the leaf
#' physiognomic space represented in the calibration dataset of Peppe et al. (2011).
#' The fossil sites being tested are placed along the CCA axes.  If a fossil site
#' falls outside of the plotted calibration space, paleoclimate reconstructions
#' for that fossil site should be treated with caution.
#'
#' @param dilp_table The results of a call to [dilp()]
#' @param physiognomy_calibration A physiognomic calibration dataset. Defaults to an internal version of
#' \code{\link{physiognomy_calibration_data}}.
#' @param climate_calibration A climate calibration dataset. Defaults to an internal version of
#' \code{\link{climate_calibration_data}}.
#'
#' @references
#' * Peppe, D.J., Royer, D.L., Cariglino, B., Oliver, S.Y., Newman, S., Leight, E., Enikolopov, G., Fernandez-Burgos, M., Herrera, F., Adams, J.M., Correa, E., Currano, E.D., Erickson, J.M., Hinojosa, L.F., Hoganson, J.W., Iglesias, A., Jaramillo, C.A., Johnson, K.R., Jordan, G.J., Kraft, N.J.B., Lovelock, E.C., Lusk, C.H., Niinemets, Ü., Peñuelas, J., Rapson, G., Wing, S.L. and Wright, I.J. (2011), Sensitivity of leaf size and shape to climate: global patterns and paleoclimatic applications. New Phytologist, 190: 724-739. https://doi.org/10.1111/j.1469-8137.2010.03615.x
#' @return A ggplot2 plot
#' @export
#'
#' @examples
#' results <- dilp(McAbeeExample)
#' dilp_cca(results)
#'
dilp_cca <- function(dilp_table, physiognomy_calibration = physiognomy_calibration_data, climate_calibration = climate_calibration_data, colorby = "data") {
  # sort by alphabetical order, to ensure sites line up between this and the climate dataframe
  physiognomy_calibration <- physiognomy_calibration[order(physiognomy_calibration$Site), ]
  colnames(physiognomy_calibration) <- colnameClean(physiognomy_calibration)
  colnames(climate_calibration) <- colnameClean(climate_calibration)

  # remove the site column
  cca_leaf <- physiognomy_calibration %>%
    dplyr::select(-"site")

  # sort by alphabetical order, to ensure sites line up between this and the leaf trait dataframe
  climate_calibration <- climate_calibration[order(climate_calibration$site), ]

  # remove the site column
  cca_climate <- climate_calibration %>%
    dplyr::select(-c("site", "koppen", "whittaker"))

  ###### Preform the CCA, excluding the fossil site(s)
  cca_analysis <- vegan::cca(X = cca_leaf, Y = cca_climate, scale = "TRUE")

  ##### Bring in the fossil sites(s) passively as an unknown sample
  #### Create new data frame that has the fossil site(s)
  coln <- colnames(physiognomy_calibration)
  cca_fossil_site <- subset(dilp_table$processed_site_data, select = coln)
  cca_fossil <- cca_fossil_site %>%
    dplyr::select(-"site")

  ## Predict CCA1 and CCA2 scores from the existing model
  cca_analysis_fossil <- stats::predict(cca_analysis, type = "wa", newdata = cca_fossil)


  ######### Plot CCA
  #### Convert CCA1 and CC2 axes scores to a dataframe
  cca_df <- data.frame(vegan::scores(cca_analysis, display = "sites"))
  # designate as calibration data
  cca_df$data <- "calibration"
  # bring the site column back
  cca_df <- cbind(cca_df, dplyr::select(climate_calibration, "site", "koppen", "whittaker"))

  #### Merge fossil data
  # Convert CCA predictions to dataframe
  cca_df_fossil <- as.data.frame(cca_analysis_fossil)
  # designate as fossil data
  cca_df_fossil$data <- "fossil"
  # bring the site column back
  cca_df_fossil <- cbind(cca_df_fossil, dplyr::select(dilp_table$processed_site_data, "site"))
  cca_df_fossil$koppen <- "Fossil"
  cca_df_fossil$whittaker <- "Fossil"
  # merge the fossil and calibration CCA data
  cca_df_all <- rbind(cca_df, cca_df_fossil)

  #### Create the plot
  # Create the convex hull
  hull_data <- cca_df %>% dplyr::slice(grDevices::chull(.data$CCA1, .data$CCA2))

  colorby <- colorby %>% stringr::str_to_lower()
  if(!(colorby %in% colnames(cca_df_all))) {
    colorby <- "data"
  }

  colnames(cca_df_all)[colnames(cca_df_all) == colorby] <- "colorby"

    # Plot
  cca_plot <- ggplot2::ggplot(data = cca_df_all, ggplot2::aes(x = .data$CCA1, y = .data$CCA2)) +
    ggplot2::geom_point(ggplot2::aes(color = .data$colorby, shape = .data$data), size = 4) +
    ggplot2::theme_classic() +
    ggplot2::geom_polygon(
      data = hull_data,
      ggplot2::aes(fill = data),
      alpha = 0.2,
      show.legend = FALSE
    ) +
    ggrepel::geom_label_repel(
      data = cca_df_fossil, ggplot2::aes(label = .data$site),
      box.padding = 0.35, point.padding = 0.5, segment.color = "grey50", max.overlaps = 50
    ) +
    scale_color_manual()
  return(cca_plot)
}
