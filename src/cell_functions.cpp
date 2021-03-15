//#define ARMA_64BIT_WORD 1  // <-- this causes compilation errors
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;


enum class FireStatus {
  UNDEFINED = -1,
  Unknown = 0,
  NoFireRegime = 1,
  TooFrequentlyBurnt = 2,
  Vulnerable = 3,
  LongUnburnt = 4,
  WithinThreshold = 5
};


// Determine fire status given interval thresholds and fire history
//
// This function determines the fire status of a location (e.g. point or
// raster cell) given a vector of fire years, and the minimum and maximum
// tolerable fire intervals (years) for the vegetation at the location.
//
// @param fireyears Vector of integer fire years. Any years beyond the
//   current query year are ignored, as are duplicate years.
//
// @param min_threshold Minimum tolerable fire interval (years).
//
// @param max_threshold Maximum tolerable fire interval (years).
//
// @param query_years One or more reference years for which to report
//   fire status.
//
// @param base_year Default start year of the data period. If no fires have
//   occurred, and \code{estimate_base_tsf == true} a time since fire value for
//   the query year is set by assuming that a fire occurred N years before the
//   base_year, where N is the mid-point of the vegetation's tolerable fire
//   interval. If fires have occurred, but only after the base year, then we
//   assume an earlier fire based on the same mid-point calculation. This only
//   applies if the minimum and maximum threshold values are not 0 (no fire
//   regime) or 9999 (vegetation that should not be burnt.)
//
// @param estimate_base_tsf If \code{true} and no fire years are earlier than
//   \code{base_year}, a pre-base year fire date is assumed based on the
//   mid-point of interval thresholds (see above). This only applies if
//   the minimum and maximum threshold values are not 0 (no fire regime) or
//   9999 (vegetation that should not be burnt.)
//
// @param quiet If true (default) do not print debugging output.
//
// @return A vector giving fire status at each query year.
//
// [[Rcpp::export]]
arma::irowvec do_cell_fire_status(const arma::ivec& fireyears,
                                  const int min_threshold,
                                  const int max_threshold,
                                  const arma::ivec& query_years,
                                  const int base_year,
                                  const bool estimate_base_tsf,
                                  const bool quiet = true) {

  enum class IntervalStatus {
    UNDEFINED = -1,
    WithinThreshold = 1,
    Vulnerable = 2,
    TooFrequent = 3
  };

  arma::irowvec rtn_status(query_years.n_elem, arma::fill::zeros);

  // --- Short-cuts for special cases
  FireStatus special_status = FireStatus::UNDEFINED;

  if (max_threshold == 0 && min_threshold == 0) {
    // Flag values for no fire regime defined
    special_status = FireStatus::NoFireRegime;

  } else if (min_threshold == 9999 && max_threshold == 9999) {
    // Flag values indicating a cell with veg that should not be burnt
    special_status = (fireyears.n_elem > 0) ? FireStatus::TooFrequentlyBurnt : FireStatus::Vulnerable;
  }

  if (special_status != FireStatus::UNDEFINED) {
    rtn_status.fill(static_cast<int>(special_status));

  } else {
    // --- General cases
    //
    for (unsigned int iyear = 0; iyear < query_years.n_elem; ++iyear) {
      int the_year = query_years(iyear);

      FireStatus fStatus = FireStatus::UNDEFINED;
      IntervalStatus iStatus = IntervalStatus::UNDEFINED;

      // Vector for query fireyears - allow space for an additional
      // year that might be added (below) before the base year
      //
      arma::uvec include_fires = arma::find(fireyears <= the_year);
      arma::ivec tmp_fireyears(include_fires.n_elem + 1);
      tmp_fireyears(0) = -1;

      if (include_fires.n_elem > 0) {
        tmp_fireyears(arma::span(1, include_fires.n_elem)) = fireyears(include_fires);
      }

      // Account for left-censoring if requested.
      // Note: this step is extra to the FireTools logic
      //
      bool leading_year = false;
      if (estimate_base_tsf &&
          (fireyears.n_elem == 0 || fireyears.min() > base_year) ) {
        int adj = (max_threshold - min_threshold) / 2;
        tmp_fireyears(0) = std::max(0, base_year - adj);
        leading_year = true;
      }

      // Subset to unique fire years in ascending order
      unsigned int ntmp = tmp_fireyears.n_elem;
      if (!leading_year) ntmp -= 1;
      arma::ivec query_fireyears = arma::unique( tmp_fireyears.tail(ntmp) );

      if (!quiet) Rcout << "fire years: " << query_fireyears << "\n";

      int FireFrequency = query_fireyears.n_elem;

      // Note: fire frequency will always be >= 1 if `estimate_base_tsf` is true
      // because a leading pre-base year will have been inserted above.
      int TSF = (FireFrequency > 0) ?
        (the_year - query_fireyears.max()) : (the_year - base_year + 1);

      // --- Decision steps begin ---
      //
      if (FireFrequency == 0) {
        fStatus = TSF > max_threshold ? FireStatus::LongUnburnt : FireStatus::Unknown;

      } else {
        // Examine fire intervals
        //
        arma::ivec intervals = arma::diff(query_fireyears);

        iStatus = IntervalStatus::WithinThreshold;

        if (!quiet) Rcout << "Assessing intervals: \n";

        for(unsigned int i_intv = 0; i_intv < intervals.n_elem; ++i_intv) {
          int the_interval = intervals(i_intv);

          if (the_interval < min_threshold) {
            if (iStatus == IntervalStatus::WithinThreshold) {
              iStatus = IntervalStatus::Vulnerable;
            } else {
              iStatus = IntervalStatus::TooFrequent;
            }
          } else if (the_interval > 2 * min_threshold) {
            iStatus = IntervalStatus::WithinThreshold;

          } else if(iStatus == IntervalStatus::WithinThreshold || iStatus == IntervalStatus::Vulnerable) {
            iStatus = IntervalStatus::WithinThreshold;
          }

          if (!quiet) Rcout << "  " << the_interval << "yr status=" << static_cast<int>(iStatus) << "\n";
        }

        if (iStatus == IntervalStatus::TooFrequent) {
          fStatus = TSF > 2 * min_threshold ?
          FireStatus::WithinThreshold : FireStatus::TooFrequentlyBurnt;

        } else if (TSF < min_threshold) {
          fStatus = FireStatus::Vulnerable;

        } else if (TSF > max_threshold) {
          fStatus = FireStatus::LongUnburnt;

        } else {
          fStatus = FireStatus::WithinThreshold;
        }
      }

      if ( fStatus == FireStatus::UNDEFINED ) {
        fStatus = FireStatus::Unknown;
      }

      rtn_status(iyear) = static_cast<int>(fStatus);
    }
  }

  return rtn_status;
}


// Determine fire status for a set of locations
//
// This function takes fire history and vegetation data for a set of locations
// and calls do_cell_fire_status for each location in turn.
//
// @param firehistory A matrix or data frame with two columns: integer
//   location ID (e.g. raster cell number); and integer fire year.
//
// @param veg A matrix or data frame with two columns: integer location ID
//   (e.g. raster cell number); and integer vegetation type code. There must
//   only be one row for each location ID.
//
// @param thresholds A matrix or data frame with three columns: integer veg
//   code; minimum tolerable interval (years); maximum tolerable interval
//   (years). Setting both the minimum and maximum values to zero indicates that
//   the vegetation type has no fire regime defined. Setting both values to
//   \code{NA} indicates that the vegetation type should never be burnt.
//
// @param query_years One or more reference years for which to report
//   fire status.
//
// @param base_year Default start year of the data period. If no fires have
//   occurred, and \code{estimate_base_tsf == true} a time since fire value for
//   the query year is set by assuming that a fire occurred N years before the
//   base_year, where N is the mid-point of the vegetation's tolerable fire
//   interval. If fires have occurred, but only after the base year, then we
//   assume an earlier fire based on the same mid-point calculation. This only
//   applies if the minimum and maximum threshold values are not 0 (no fire
//   regime) or 9999 (vegetation that should not be burnt.)
//
// @param estimate_base_tsf If \code{true} and no fire years are earlier than
//   \code{base_year}, a pre-base year fire date is assumed based on the
//   mid-point of interval thresholds (see above). This only applies if
//   the minimum and maximum threshold values are not 0 (no fire regime) or
//   9999 (vegetation that should not be burnt.)
//
// @param quiet If true (default) do not print debugging output.
//
// @return
// [[Rcpp::export]]
arma::imat table_fire_status(const arma::imat& firehistory,
                             const arma::imat& veg,
                             const arma::imat& veg_thresholds,
                             const arma::ivec& query_years,
                             const int base_year,
                             const bool estimate_base_tsf,
                             const bool quiet = true) {

  arma::ivec cells = arma::unique( veg.col(0) );

  // output matrix has cell column, then a column for each query year
  const unsigned int Nyears = query_years.n_elem;
  arma::imat fire_status(cells.n_elem, 1+Nyears);

  for (unsigned int i = 0; i < cells.n_elem; ++i) {
    // veg unit
    arma::uvec rr = arma::find(veg.col(0) == cells(i));
    int iveg = veg(rr(0), 1);

    // thresholds
    rr = arma::find(veg_thresholds.col(0) == iveg);
    int min_threshold = veg_thresholds(rr(0), 1);
    int max_threshold = veg_thresholds(rr(0), 2);

    // fire years
    rr = arma::find(firehistory.col(0) == cells(i));
    arma::ivec fireyears(rr.n_elem);
    for (unsigned int k = 0; k < rr.n_elem; ++k) {
      fireyears(k) = firehistory(rr(k), 1);
    }

    fire_status(i, 0) = cells(i);

    fire_status(i, arma::span(1,Nyears)) = do_cell_fire_status(
                fireyears,
                min_threshold, max_threshold,
                query_years, base_year, quiet);
  }

  return fire_status;
}

