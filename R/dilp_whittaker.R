base_whittaker <- function(){
  biomes <- Whittaker_biomes
  Whittaker_biomes$biomeorder <- 1:nrow(Whittaker_biomes)
  biome<- c("Tundra",
                    "Boreal forest",
                    "Temperate seasonal forest",
                    "Temperate rain forest",
                    "Tropical rain forest",
                    "Tropical seasonal forest/savanna",
                    "Subtropical desert",
                    "Temperate grassland/desert",
                    "Woodland/shrubland")
  order <- c("1", "2", "3", "4", "5", "6", "7", "8", "9")

  colors <- c("#C1E1DD", "#A5C790", "#97B669", "#75A95E", "#317A22", "#A09700", "#DCBB50", "#FCD57A", "#D16E3F")

  all_biomes <- data.frame(biome, order)
  all_biomes <- merge(all_biomes, Whittaker_biomes)
  all_biomes <- all_biomes[order(all_biomes$order, all_biomes$biomeorder),]
  all_biomes$biome <- factor(all_biomes$biome, levels = biome)
  plot <- ggplot2::ggplot() + ggplot2::labs(x = paste0("Mean Annual Temperature (", "\u00B0", "C)"),
                                            y = "Mean Annual Precipitation (cm/yr)") +
    ggplot2::geom_polygon(data = all_biomes, ggplot2::aes(x = .data$temp_c, y = .data$precp_cm, fill = biome),
                          color = "black", linewidth = 0.5) +
    ggplot2::scale_fill_manual(values = colors, name = "Biome") +
    ggplot2::theme_bw()
  return(plot)
}

#' Plot climate reconstructions on a Whittaker Biome plot
#' @description
#' `dilp_whittaker()` plots `dilp()` outputs onto a Whittaker Biome plot.  Base
#' Whittaker Plot from the [plotbiomes](https://github.com/valentinitnelav/plotbiomes)
#' package by Ștefan Valentin and Sam Levin.
#'
#' @param climate_data A data frame containing either the direct output of a
#' [dilp()] or [dilp_pgls()] call, or the `$results` element from either.
#' Can also be a data frame with the following required columns (MLR or PIP):
#' * `site`
#' * `MAT.MLR` or `MAT.PIP`
#' * `MAT.MLR.error` or `MAT.PIP.error`
#' * `MAP.MLR` or `MAP.PIP`
#' * `MAP.MLR.error.minus` / `MAP.MLR.error.plus`, or
#'   `MAP.PIP.error.minus` / `MAP.PIP.error.plus`
#'
#'
#' @returns A modifiable ggplot with dilp climate-reconstructed sites plotted
#' onto a Whittaker diagram.
#'
#' @references
#' Valentin Ștefan, & Sam Levin. (2018). plotbiomes: R package for plotting Whittaker biomes with ggplot2 (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.7145245
#'
#' @export
#'
#' @examples
#' results <- dilp(McAbeeExample)
#' dilp_whittaker(results)
#'
dilp_whittaker <- function(climate_data){
  if (is.list(climate_data)) {
    if(!is.null(climate_data$results)) {
      climate_data <- climate_data$results
    }
  }

  # Detect whether data comes from dilp_pgls() (PIP columns) or dilp() (MLR columns)
  if ("MAT.PIP" %in% colnames(climate_data)) {
    mat_col       <- "MAT.PIP"
    map_col       <- "MAP.PIP"
    mat_err       <- "MAT.PIP.error"
    map_err_plus  <- "MAP.PIP.error.plus"
    map_err_minus <- "MAP.PIP.error.minus"
  } else {
    mat_col       <- "MAT.MLR"
    map_col       <- "MAP.MLR"
    mat_err       <- "MAT.MLR.error"
    map_err_plus  <- "MAP.MLR.error.plus"
    map_err_minus <- "MAP.MLR.error.minus"
  }

  base_plot <- base_whittaker()

  final_plot <- base_plot +
    ggrepel::geom_label_repel(data = climate_data,
                              ggplot2::aes(x = .data[[mat_col]], y = .data[[map_col]], label = .data$site),
                              box.padding = 0.35, point.padding = 0.5, segment.color = "grey50", max.overlaps = 50
    ) +
    ggplot2::geom_point(data = climate_data,
                        ggplot2::aes(x = .data[[mat_col]], y = .data[[map_col]]),
                        size = 3) +
    ggplot2::geom_errorbar(data = climate_data,
                           ggplot2::aes(xmin = .data[[mat_col]] - .data[[mat_err]],
                                        xmax = .data[[mat_col]] + .data[[mat_err]],
                                        y = .data[[map_col]])) +
    ggplot2::geom_errorbar(data = climate_data,
                           ggplot2::aes(ymin = .data[[map_col]] - .data[[map_err_minus]],
                                        ymax = .data[[map_col]] + .data[[map_err_plus]],
                                        x = .data[[mat_col]]))

  return(final_plot)
}
