#' Fire interval threshold analysis for one or more point locations or raster cells
#'
#' @param firehistory A matrix or data frame with two columns: integer
#'   location ID (e.g. raster cell number); and integer fire year.
#'
#' @param veg A matrix or data frame with two columns: integer location ID
#'   (e.g. raster cell number); and integer vegetation type code. There must
#'   only be one row for each location ID.
#'
#' @param thresholds A matrix or data frame with three columns: integer veg
#'   code; minimum tolerable interval (years); maximum tolerable interval
#'   (years). Setting both the minimum and maximum values to zero indicates that
#'   the vegetation type has no fire regime defined. Setting both values to
#'   \code{NA} indicates that the vegetation type should never be burnt.
#'
#' @param query_years One or more reference years for which to report fire
#'   status.
#'
#' @param base_year The earliest year covered by the provided fire history. This
#'   can be earlier than any actual fire years specified for locations.
#'
#' @export
#'
get_fire_status <- function(firehistory,
                            veg,
                            thresholds,
                            query_years, base_year) {

  veg <- as.matrix(veg)
  nv <- nrow(veg)

  if (nv == 0) {
    # Nothing to do - return an empty matrix
    out <- matrix(0, nrow = 0, ncol = 2, dimnames = list(NULL, c("cell", "status")))
    return(out)
  }

  if (anyNA(veg)) {
    stop("one or more missing values in veg table")
  }

  if (nv > length(unique(veg[,1]))) {
    stop("veg table has one or more repeated location IDs")
  }

  firehistory <- as.matrix(firehistory)
  if (!all(firehistory[,1] %in% veg[,1])) {
    stop("one or more location IDs in firehistory table are not in the veg table")
  }

  # Remove any missing fire years
  keep <- !is.na(firehistory[,2])
  firehistory <- firehistory[keep, , drop=FALSE]

  thresholds <- as.matrix(thresholds)

  # Replace NA thresholds with flag value 9999
  namin <- is.na(thresholds[,2])
  namax <- is.na(thresholds[,3])
  if (!all(namin == namax)) {
    stop("One or more records have NA for one but not both thresholds")
  }
  thresholds[namin, 2] <- 9999
  thresholds[namax, 3] <- 9999

  # Check no threshold pairs are in the wrong order
  if (!all(veg[,2] %in% thresholds[,1])) {
    stop("One or more veg codes in the veg table are not in the thresholds table")
  }

  if (any(thresholds[,2] > thresholds[,3])) {
    stop("One or more minimum threshold values are greater than the corresponding maximum value")
  }

  if (!all(query_years > base_year)) {
    stop("query years must be later than the specified base year")
  }

  # Call Rcpp function and return the resulting matrix
  table_fire_status(firehistory, veg, thresholds, query_years, base_year)
}
