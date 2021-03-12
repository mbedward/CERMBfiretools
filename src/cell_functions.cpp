#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

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


// Helper function to find the first row index where a given column has
// the search value.
//
// [[Rcpp::export]]
int find_row(IntegerMatrix x, int col_index, int value) {
  IntegerVector c = x.column(col_index);
  int *i = std::find(c.begin(), c.end(), value);
  if (i == c.end()) return -1;
  else return i - c.begin();
}

// Helper function to get the subset of a matrix where column col_index
// is equal to value.
//
// [[Rcpp::export]]
IntegerMatrix mat_subset(IntegerMatrix x, int col_index, int value) {
  IntegerMatrix::Column c = x.column(col_index);
  LogicalVector b = (c == value);

  int nr = sum(b);
  int nc = x.ncol();

  IntegerMatrix res(nr, nc);

  for (int ir = 0, ib = 0; ib < x.nrow(); ++ib) {
    if (b[ib] == TRUE) {
      res(ir++, _) = x(ib, _);
    }
  }

  return res;
}


//' Determine fire status given interval thresholds and fire history
//'
//' This function determines the fire status of a location (e.g. point or
//' raster cell) given a vector of fire years, and the minimum and maximum
//' tolerable fire intervals (years) for the vegetation at the location.
//'
//' @param fireyears Vector of integer fire years. Any years beyond year_query
//'   are ignored, as are duplicate years.
//'
//' @param min_threshold Minimum tolerable fire interval (years).
//'
//' @param max_threshold Maximum tolerable fire interval (years).
//'
//' @param year_query Reference year for the analysis.
//'
//' @param year_base Year taken as the start of the data period.
//'
//' @param quiet If true (default) do not print debugging output.
//'
//' @return An integer code for fire status.
//'
//'
//' @export
// [[Rcpp::export]]
int cell_fire_status(IntegerVector fireyears,
                     int min_threshold, int max_threshold,
                     int year_query, int year_base,
                     bool quiet = true) {

  enum class IntervalStatus {
    UNDEFINED = -1,
    WithinThreshold = 1,
    Vulnerable = 2,
    TooFrequent = 3
  };

  FireStatus fStatus = FireStatus::UNDEFINED;
  IntervalStatus iStatus = IntervalStatus::UNDEFINED;

  // Remove any NA values and years beyond year_query from the
  // fireyears vector
  fireyears = ifelse(fireyears > year_query, NA_INTEGER, fireyears);
  fireyears = na_omit(fireyears);

  // Ensure years are sorted and remove any duplicates
  fireyears = unique(fireyears.sort());

  // Account for left-censoring (this step is extra to the FireTools logic)
  if (fireyears.length() == 0 || fireyears(0) > year_base) {
    int adj = (max_threshold - min_threshold) / 2;
    fireyears.push_front(year_base - adj);
  }

  if (!quiet) Rcout << "fire years: " << fireyears << "\n";


  int FireFrequency = fireyears.length();
  int TSF = FireFrequency > 0 ? (year_query - max(fireyears)) : (year_query - year_base + 1);

  if (IntegerVector::is_na(min_threshold) && IntegerVector::is_na(max_threshold)) {
    // Cell with veg that should not be burnt
    fStatus = FireFrequency > 0 ? FireStatus::TooFrequentlyBurnt : FireStatus::Vulnerable;

  } else if(max_threshold == 0 && min_threshold == 0) {
    fStatus = FireStatus::NoFireRegime;

  } else if (FireFrequency == 0) {
    fStatus = TSF > max_threshold ? FireStatus::LongUnburnt : FireStatus::Unknown;

  } else {
    // None of the special cases above applied, so we do the
    // interval analysis
    IntegerVector intervals = diff(fireyears);

    iStatus = IntervalStatus::WithinThreshold;

    if (!quiet) Rcout << "Assessing intervals: \n";
    for(IntegerVector::iterator intv = intervals.begin(); intv != intervals.end(); ++intv) {
      if (*intv < min_threshold) {
        if (iStatus == IntervalStatus::WithinThreshold) {
          iStatus = IntervalStatus::Vulnerable;
        } else {
          iStatus = IntervalStatus::TooFrequent;
        }
      } else if (*intv > 2 * min_threshold) {
        iStatus = IntervalStatus::WithinThreshold;

      } else if(iStatus == IntervalStatus::WithinThreshold || iStatus == IntervalStatus::Vulnerable) {
        iStatus = IntervalStatus::WithinThreshold;
      }

      if (!quiet) Rcout << "  " << *intv << "yr status=" << static_cast<int>(iStatus) << "\n";
    }

    if (iStatus == IntervalStatus::TooFrequent) {
      fStatus = TSF > 2 * min_threshold ? FireStatus::WithinThreshold : FireStatus::TooFrequentlyBurnt;

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

  return static_cast<int>(fStatus);
}


//' Determine fire status for a set of locations
//'
//' @export
// [[Rcpp::export]]
IntegerMatrix table_fire_status(IntegerMatrix cell_firehistory,
                                IntegerMatrix cell_veg,
                                IntegerMatrix veg_thresholds,
                                int year_query, int year_base,
                                bool quiet = true) {

  IntegerVector cells = unique( cell_veg(_, 0) ).sort();
  if (cells.length() == 0) return IntegerMatrix(0, 2); // nothing to do

  IntegerMatrix status( cells.length(), 2 );

  int k = 0;
  for (IntegerVector::iterator icell = cells.begin(); icell != cells.end(); ++icell, ++k) {
    int vr = find_row(cell_veg, 0, *icell);
    int veg = cell_veg(vr, 1);

    int tr = find_row(veg_thresholds, 0, veg);
    int min_threshold = veg_thresholds(tr, 1);
    int max_threshold = veg_thresholds(tr, 2);

    IntegerMatrix fmat = mat_subset(cell_firehistory, 0, *icell);
    IntegerVector fireyears = fmat(_, 1);

    status(k, 0) = *icell;
    status(k, 1) = cell_fire_status(fireyears, min_threshold, max_threshold,
                                    year_query, year_base, quiet);
  }

  return status;
}

