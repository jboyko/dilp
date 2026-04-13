#' Estimate temperature with simple linear regression
#'
#' @description
#' `temp_slr()` will produce estimates of mean annual temperature and standard error
#' using leaf margin analysis. There are different ways to represent error. The most simple is using the standard error of the regression.
#' These are listed in the table below. However, this is not the only source of uncertainty and is too simplistic a measure of error.
#' This function instead uses the method outlined in Miller et al. 2006, and reported in Peppe et al 2018 (eq. 4), which also accounts
#' for binomial sampling error and overdispersion, offering what we consider a best practice approach. The standard error of the regression provides a minimum error value.
#' Note, Peppe et al. 2018 suggests that a conservative minimum uncertainty for all leaf margin analysis results is probably +/- 5 degrees Celsius.
#'
#' Standard error of regression:
#'
#' | Regression      | SE |
#' | --------------- | --- |
#' | Peppe2018       | 4.5 |
#' | Peppe2011       | 4.8 |
#' | Peppe2011NH     | 3.4 |
#' | Miller2006      |  -  |
#' | WingGreenwood   | 0.8 |
#' | Wilf1997        | 2.0 |
#' | KowalskiDilcher | 3.6 |
#' @md
#'
#' @param data A data frame that must include the columns "morphotype" and "margin".
#' Can be leaf or species level data.
#' @param regression A string representing one of the following pre-loaded regressions:
#' * "Peppe2018" - for global temperature estimates
#' * "Peppe2011" - The Americas, Japan, and Oceania
#' * "Peppe2011NH" - Peppe 2011 (Northern Hemisphere only)
#' * "Miller2006" - North and Central America
#' * "WingGreenwood" - East Asia - original leaf margin analysis regression
#' * "Wilf1997" - The Americas
#' * "KowalskiDilcher" - North America
#'
#' @param slope Slope, if using a custom regression
#' @param constant Constant, if using a custom regression
#' @param error Standard error, if using a custom regression
#'
#' @return A table with MAT estimates for each site
#' @references
#' * Kowalski, E.A. & Dilcher, D.L. (2003). Warmer paleotemperatures for terrestrial ecosystems. Proceedings of the National Academy of Sciences, 100, 167–170.
#' * Miller, I. M., Brandon, M. T., & Hickey, L. J. (2006). Using leaf margin analysis to estimate mid-Cretaceous (Albian) paleolatitude of the Baja BC block. Earth and Planetary Science Letters, 245, 95–114.
#' * Peppe, D. J., Baumgartner, A., Flynn, A., & Blonder, B. (2018). Reconstructing paleoclimate and paleoecology using fossil leaves. Methods in paleoecology: Reconstructing Cenozoic terrestrial environments and ecological communities, 289-317.
#' * Peppe, D.J., Royer, D.L., Cariglino, B., Oliver, S.Y., Newman, S., Leight, E., Enikolopov, G., Fernandez-Burgos, M., Herrera, F., Adams, J.M., Correa, E., Currano, E.D., Erickson, J.M., Hinojosa, L.F., Hoganson, J.W., Iglesias, A., Jaramillo, C.A., Johnson, K.R., Jordan, G.J., Kraft, N.J.B., Lovelock, E.C., Lusk, C.H., Niinemets, Ü., Peñuelas, J., Rapson, G., Wing, S.L. and Wright, I.J. (2011), Sensitivity of leaf size and shape to climate: global patterns and paleoclimatic applications. New Phytologist, 190: 724-739. https://doi.org/10.1111/j.1469-8137.2010.03615.x
#' * Wing, S., & Greenwood, D. R. (1993). Fossils and fossil climate: the case for equable continental interiors in the Eocene. Philosophical Transactions of the Royal Society of London Series B, 341, 243–252.
#' * Wilf, P. (1997). When are leaves good thermometers? A new case for leaf margin analysis. Paleobiology, 23, 373–390.
#'
#' @export
#'
#' @examples
#' temp_slr(McAbeeExample, regression = "Peppe2011")
temp_slr <- function(data, regression = "Peppe2018", slope = NULL, constant = NULL, error = NULL) {
  cleaned <- colnameClean(data)
  if (!"morphotype" %in% cleaned) {
    cleaned <- stringr::str_replace_all(cleaned, "^species$", "morphotype")
  }
  colnames(data) <- cleaned
  required_columns <- c("morphotype", "margin")
  missing_columns <- required_columns[!required_columns %in% colnames(data)]
  if (length(missing_columns) > 0) {
    stop(paste("data is missing required columns:", stringr::str_flatten(missing_columns, collapse = ", ")))
  }
  if (length(unique(data$margin)) > 2) {
    data$margin[data$margin > 0 & data$margin < 1] <- 0.5
    if (length(unique(data$margin)) > 2) {
      warning("Margin states outside the bounds of [0 - 1] present")
    }
  }
  if (is.null(c(slope, constant, error))) {
    regression <- grab_regression(regression, "temp")
    slope <- regression$slope
    constant <- regression$constant
    error <- regression$error
  } else if (is.null(slope) | !is.numeric(slope)) {
    stop("Attempting to use custom regression but slope not a valid number")
  } else if (is.null(constant) | !is.numeric(constant)) {
    stop("Attempting to use custom regression but constant not a valid number")
  } else if (is.null(error) | !is.numeric(error)) {
    stop("Attempting to use custom regression but error not a valid number")
  }

  sites <- data %>%
    dplyr::distinct(.data$site) %>%
    cbind("n" = NA, "lower" = NA, "MAT" = NA, "upper" = NA)

  for (i in 1:nrow(sites)) {
    site_subset <- dplyr::filter(data, data$site == as.character(sites$site[i])) %>%
      dplyr::group_by(.data$morphotype) %>%
      dplyr::summarize(margin = mean(.data$margin, na.rm = TRUE)) %>%
      stats::na.omit()

    num_morph <- nrow(site_subset)
    margin_p <- (100 * sum(site_subset$margin) / num_morph)
    miller_error <- miller_error(slope, num_morph, margin_p)

    if(miller_error > error){
      error <- miller_error
    }

    value <- slope * margin_p + constant
    sites$n[i] <- num_morph
    sites$MAT[i] <- value
    sites$lower[i] <- value - error
    sites$upper[i] <- value + error
  }
  return(sites)
}

miller_error <- function(slope, num_morph, margin_p) {
  margin_p <- margin_p / 100
  error <- (100*slope)*(sqrt((1+0.052*(num_morph-1)*margin_p*(1-margin_p))*((margin_p*(1-margin_p))/num_morph)))
  return(error)
}


#' Estimate precipitation with simple linear regression
#'
#' @description
#' `precip_slr()` will produce estimates of mean annual precipitation and standard error
#' using leaf area analysis.
#' @param data A data frame that must include the columns "morphotype", "leaf_area", and "specimen_number".
#' Must be leaf level data.
#' @param regression A string representing one of the following pre-loaded regressions:
#' * "Peppe2018" - for global precipitation estimates
#' * "Peppe2011" - The Americas, Japan, and Oceania
#' * "Jacobs2002" - Africa
#' * "Wilf1998" - The Americas and Africa
#'
#' @param slope Slope, if using a custom regression
#' @param constant Constant, if using a custom regression
#' @param error Standard error, if using a custom regression
#'
#' @return A table with MAP estimates for each site
#' @references
#' * Peppe, D. J., Baumgartner, A., Flynn, A., & Blonder, B. (2018). Reconstructing paleoclimate and paleoecology using fossil leaves. Methods in paleoecology: Reconstructing Cenozoic terrestrial environments and ecological communities, 289-317.
#' * Peppe, D.J., Royer, D.L., Cariglino, B., Oliver, S.Y., Newman, S., Leight, E., Enikolopov, G., Fernandez-Burgos, M., Herrera, F., Adams, J.M., Correa, E., Currano, E.D., Erickson, J.M., Hinojosa, L.F., Hoganson, J.W., Iglesias, A., Jaramillo, C.A., Johnson, K.R., Jordan, G.J., Kraft, N.J.B., Lovelock, E.C., Lusk, C.H., Niinemets, Ü., Peñuelas, J., Rapson, G., Wing, S.L. and Wright, I.J. (2011), Sensitivity of leaf size and shape to climate: global patterns and paleoclimatic applications. New Phytologist, 190: 724-739. https://doi.org/10.1111/j.1469-8137.2010.03615.x
#' * Jacobs, B. F. (2002). Estimation of low-latitude paleoclimates using fossil angiosperm leaves: examples from the Miocene Tugen Hills, Kenya. Paleobiology, 28, 399–421.
#' * Wilf, P. (2008). Fossil angiosperm leaves: paleobotany’s difficult children prove themselves. Paleontological Society Papers, 14, 319–333.
#' @export
#'
#' @examples
#' precip_slr(McAbeeExample, regression = "Peppe2011")
precip_slr <- function(data, regression = "Peppe2018", slope = NULL, constant = NULL, error = NULL) {
  colnames(data) <- colnameClean(data)

  if ("leaf_area" %in% colnames(data)) {

  } else if ("blade_area" %in% colnames(data)) {
    if ("petiole_area" %in% colnames(data)) {

    } else {
      data[["petiole_area"]] <- NA
      warning("petiole_area missing, blade_area used as leaf_area", call. = FALSE)
    }
    data$petiole_area <- ifelse(data$blade_area > 0 & is.na(data$petiole_area), 0, data$petiole_area)
    data$leaf_area <- data$blade_area + data$petiole_area
  }
  required_columns <- c("morphotype", "leaf_area", "specimen_number")
  missing_columns <- required_columns[!required_columns %in% colnames(data)]
  if (length(missing_columns) > 0) {
    stop(paste("data is missing required columns:", stringr::str_flatten(missing_columns, collapse = ", ")))
  }
  if (is.null(c(slope, constant, error))) {
    regression <- grab_regression(regression, "precip")
    slope <- regression$slope
    constant <- regression$constant
    error <- regression$error
  } else if (is.null(slope) | !is.numeric(slope)) {
    stop("Attempting to use custom regression but slope not a valid number")
  } else if (is.null(constant) | !is.numeric(constant)) {
    stop("Attempting to use custom regression but constant not a valid number")
  } else if (is.null(error) | !is.numeric(error)) {
    stop("Attempting to use custom regression but error not a valid number")
  }

  sites <- data %>%
    dplyr::distinct(.data$site) %>%
    cbind("n" = NA, "lower" = NA, "MAP" = NA, "upper" = NA)

  for (i in 1:nrow(sites)) {
    site_subset <- dplyr::filter(data, data$site == as.character(sites$site[i])) %>%
      dplyr::group_by(.data$morphotype) %>%
      dplyr::mutate(ln_leaf_area = log(100 * .data$leaf_area)) %>%
      dplyr::summarize(ln_leaf_area = mean(.data$ln_leaf_area, na.rm = TRUE)) %>%
      stats::na.omit()

    value <- (mean(site_subset$ln_leaf_area) * slope) + constant
    sites$n[i] <- nrow(site_subset)
    sites$MAP[i] <- exp(value)
    sites$lower[i] <- exp(value) - exp(value - error)
    sites$upper[i] <- exp(value) + exp(value + error)
  }
  return(sites)
}
