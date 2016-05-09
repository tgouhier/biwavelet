#include <stdlib.h>
#include <R_ext/Arith.h>
#include "pivot.h"
#include "quantile.h"

/**
 * Here I first put all NAs to the end, then call the pivot function to find
 * the appropriate quantile of the remaining (finite) entries.
 *
 * q is the quantile: 1/2 will give exactly the median above.
 *
 * General notes about handling missing data, zero MAD etc:
 * The idea is that bicor should switch to cor whenever it is feasible,
 * it helps, and it is requested:
 * (1) if median is NA, the mean would be NA as well, so there's no point in
 *     switching to Pearson
 * (2) In the results, columns and rows corresponding to input with NA
 *     means/medians are NA'd out.
 * (3) The convention is that if zeroMAD is set to non-zero, it is the index of
 *     the column in which MAD is zero (plus one for C indexing)
 */
double quantile(double *x, const size_t n, const double q) {

  // Put all NA's at the end.
  size_t bound = n;
  for (size_t i = n; i>0;) {
    i--;
    if (ISNAN(x[i])) {
       bound--;
       x[i] = x[bound];
       x[bound] = NA_REAL;
    }
  }

  // Any non-NA's left?
  if(bound == 0) {
    return NA_REAL;
  }

  return pivot(x, bound, 1.0 * (bound - 1) * q);
}
