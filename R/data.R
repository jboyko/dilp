#' McAbee Example Data
#'
#' Leaf physiognomic data of specimens collected from the McAbee Fossil Beds
#' in British Columbia, Canada (Lowe et al. 2018).
#'
#' @format ## `McAbeeExample`
#' A data frame with 192 rows and 18 columns:
#' \describe{
#'   \item{Site}{Stratigraphic layer or locality}
#'   \item{Specimen Number}{Repository number for individual specimen}
#'   \item{Morphotype}{Morphotype the specimen belongs to}
#'   \item{Measurer comments}{Additional notes about the specimen or its measurements}
#'   \item{Margin}{Whether the margin is toothed (0) or entire (1)}
#'   \item{Petiole Width}{The width of the petiole at the basalmost point of insertion into the leaf lamina}
#'   \item{Blade area}{The reconstructed area of the leaf lamina, not including the petiole}
#'   \item{Blade perimeter}{The length of the perimeter of the leaf lamina, not including the petiole}
#'   \item{Feret}{The diameter of a circle with the same area as the leaf lamina, not including the petiole}
#'   \item{Minimum Feret}{The longest line that can be drawn between two points on the perimeter of a selection that is perpendicular to Feret length.  Approximates blade width.}
#'   \item{Raw blade area}{The area of a leaf prepared for tooth measurements that still has its teeth.}
#'   \item{Raw blade perimeter}{The perimeter of a leaf prepared for tooth measurements that still has its teeth.}
#'   \item{Internal raw blade area}{The area of a leaf prepared for tooth measurements with teeth digitally removed.}
#'   \item{Internal raw blade perimeter}{The perimeter of a leaf prepared for tooth measurements with teeth digitally removed.}
#'   \item{Length of cut perimeter}{The total length of all segments of leaf removed from the leaf blade while removing damage during preparation of the leaf.}
#'   \item{no. of primary teeth}{The number of primary teeth along the undamaged perimeter}
#'   \item{no. of secondary teeth}{The number of secondary teeth along the undamaged perimeter}
#'
#' }
#' @references
#' * Lowe, A. J., D. R. Greenwood, C. K. West, J. M. Galloway, M. Sudermann, and T. Reichgelt. 2018. Plant community ecology and climate on an upland volcanic landscape during the Early Eocene Climatic Optimum: McAbee Fossil Beds, British Columbia, Canada. Palaeogeography, Palaeoclimatology, Palaeoecology 511: 433–448.
#' @source Lowe et al. 2018
"McAbeeExample"

#' Climate Calibration Data
#'
#' Temperature and precipitation data associated with the modern localities used to calibrate the DiLP model
#'
#' @format ## `climate_calibration_data`
#' A data frame with 92 rows and 5 columns:
#' \describe{
#'   \item{Site}{Locality name}
#'   \item{MAT}{Mean Annual Temperature (celsius)}
#'   \item{MAP}{Mean Annual Precipitation (mm)}
#'   \item{Koppen}{Koppen climate classifcations}
#'   \item{Whittaker}{Whittaker climate classifications}
#' }
#' @source Peppe et al. 2011
#' @references
#' * Peppe, D.J., Royer, D.L., Cariglino, B., Oliver, S.Y., Newman, S., Leight, E., Enikolopov, G., Fernandez-Burgos, M., Herrera, F., Adams, J.M., Correa, E., Currano, E.D., Erickson, J.M., Hinojosa, L.F., Hoganson, J.W., Iglesias, A., Jaramillo, C.A., Johnson, K.R., Jordan, G.J., Kraft, N.J.B., Lovelock, E.C., Lusk, C.H., Niinemets, Ü., Peñuelas, J., Rapson, G., Wing, S.L. and Wright, I.J. (2011), Sensitivity of leaf size and shape to climate: global patterns and paleoclimatic applications. New Phytologist, 190: 724-739. https://doi.org/10.1111/j.1469-8137.2010.03615.x
"climate_calibration_data"

#' Physiognomy Calibration Data
#'
#' Leaf physiognomic data taken from the modern localities used to calibrate the DiLP model
#'
#' @format ## `physiognomy_calibration_data`
#' A data frame with 92 rows and 12 columns:
#' \describe{
#'   \item{Site}{Locality name}
#'   \item{Leaf.area}{Average leaf area at site}
#'   \item{FDR}{Feret diameter:Feret length. Describes leaf linearity compared to a circle}
#'   \item{Perimeter.ratio}{Ratio - Raw blade perimeter:Internal raw blade perimeter}
#'   \item{TC.P}{Ratio - Tooth count:Perimeter}
#'   \item{TC.IP}{Ratio - Tooth count:Internal perimeter}
#'   \item{Avg.TA}{Average area of a primary tooth}
#'   \item{TA.BA}{Ratio - Tooth area:Blade area}
#'   \item{TA.P}{Ratio - Tooth area:Perimeter}
#'   \item{TA.IP}{Ratio - Tooth area:Internal perimeter}
#'   \item{TC.BA}{Ratio - Tooth count:Blade area}
#'   \item{Margin}{Percentage of untoothed species at the site}
#' }
#' @source Peppe et al. 2011
#' @references
#' * Peppe, D.J., Royer, D.L., Cariglino, B., Oliver, S.Y., Newman, S., Leight, E., Enikolopov, G., Fernandez-Burgos, M., Herrera, F., Adams, J.M., Correa, E., Currano, E.D., Erickson, J.M., Hinojosa, L.F., Hoganson, J.W., Iglesias, A., Jaramillo, C.A., Johnson, K.R., Jordan, G.J., Kraft, N.J.B., Lovelock, E.C., Lusk, C.H., Niinemets, Ü., Peñuelas, J., Rapson, G., Wing, S.L. and Wright, I.J. (2011), Sensitivity of leaf size and shape to climate: global patterns and paleoclimatic applications. New Phytologist, 190: 724-739. https://doi.org/10.1111/j.1469-8137.2010.03615.x
"physiognomy_calibration_data"

lma_regressions <- list(
  royer_species_mean_ma = list(
    stat = "mean",
    regression_slope = 0.382,
    y_intercept = 3.070,
    unexplained_mean_square = 0.032237,
    sample_size_calibration = 667,
    mean_log_petiole_metric_calibration = -3.011,
    sum_of_squares_calibration = 182.1,
    critical_value = 1.964
  ),
  royer_site_mean_ma = list(
    stat = "mean",
    regression_slope = 0.429,
    y_intercept = 3.214,
    unexplained_mean_square = 0.005285,
    sample_size_calibration = 25,
    mean_log_petiole_metric_calibration = -2.857,
    sum_of_squares_calibration = 5.331,
    critical_value = 2.069
  ),
  lowe_site_mean_ma = list(
    stat = "mean",
    regression_slope = 0.345,
    y_intercept = 2.954,
    unexplained_mean_square = 0.01212861,
    sample_size_calibration = 70,
    mean_log_petiole_metric_calibration = -2.902972,
    sum_of_squares_calibration = 1.154691,
    critical_value = 1.995469
  ),
  lowe_site_variance_ma = list(
    stat = "variance",
    regression_slope = 0.302,
    y_intercept = 5.028,
    unexplained_mean_square = 0.1713672,
    sample_size_calibration = 70,
    mean_log_petiole_metric_calibration = -5.97104,
    sum_of_squares_calibration = 5.085184,
    critical_value = 1.995469
  )
)

dilp_parameters <- list(
  PeppeGlobal = list(
    MAT.MLR.M = 0.21,
    MAT.MLR.FDR = 42.296,
    MAT.MLR.TC.IP = -2.609,
    MAT.MLR.constant = -16.004,
    MAT.MLR.error = 4,
    MAT.SLR.M = 0.204,
    MAT.SLR.constant = 4.6,
    MAT.SLR.error = 4.8,
    MAP.MLR.LA = 0.298,
    MAP.MLR.TC.IP = 0.279,
    MAP.MLR.PR = -2.717,
    MAP.MLR.constant = 3.033,
    MAP.MLR.SE = 0.6,
    MAP.SLR.LA = 0.283,
    MAP.SLR.constant = 2.92,
    MAP.SLR.SE = 0.61
  ),
  PeppeNH = list(
    MAT.MLR.M = 0.233,
    MAT.MLR.FDR = 0,
    MAT.MLR.TC.IP = -1.547,
    MAT.MLR.constant = 8.161,
    MAT.MLR.error = 4,
    MAT.SLR.M = 0.262,
    MAT.SLR.constant = 3.167,
    MAT.SLR.error = 2.0,
    MAP.MLR.LA = NA,
    MAP.MLR.TC.IP = NA,
    MAP.MLR.PR = NA,
    MAP.MLR.constant = NA,
    MAP.MLR.SE = NA,
    MAP.SLR.LA = NA,
    MAP.SLR.constant = NA,
    MAP.SLR.SE = NA
  )
)

temp_regressions <- list(
  Peppe2018 = list(slope = 0.194, constant = 5.884, error = 4.54), # 4.54
  Peppe2011NH = list(slope = 0.262, constant = 3.167, error = 3.4), # 2.0
  Peppe2011 = list(slope = 0.204, constant = 4.6, error = 4.8), # 4.8
  WingGreenwood = list(slope = 0.306, constant = 1.141, error = .8), # 0.80
  Wilf1997 = list(slope = 0.286, constant = 2.24, error = 2.0),
  Miller2006 = list(slope = 0.290, constant = 1.320, error = 2),
  KowalskiDilcher = list(slope = 0.363, constant = 2.223, error = 3.6)
)

precip_regressions <- list(
  Peppe2018 = list(slope = 0.346, constant = 2.404, error = 0.93),
  Peppe2011 = list(slope = 0.283, constant = 2.92, error = 0.61),
  Jacobs2002 = list(slope = 0.321, constant = 2.476, error = 0.24),
  Wilf1998 = list(slope = 0.546, constant = 0.786, error = 0.36)
)

#' Whittaker Biome Data
#'
#' Delineations of Whittaker biomes from github.com/valentinitnelav/plotbiomes
#'
#' @format ## `Whittaker_biomes`
#' A data frame with points mapping out Whittaker biome delineations.
#' @source github.com/valentinitnelav/plotbiomes
#' @references
#' *Valentin Ștefan, & Sam Levin. (2018). plotbiomes: R package for plotting Whittaker biomes with ggplot2 (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.7145245
"Whittaker_biomes"

utils::globalVariables(c(
  "climate_calibration_data",
  "physiognomy_calibration_data",
  "pgls_beta_mat",
  "pgls_beta_map",
  "pgls_lambda_mat",
  "pgls_lambda_map",
  "pgls_Ke_mat",
  "pgls_Ke_map",
  "pgls_pred_names",
  "pgls_name_table"
))
