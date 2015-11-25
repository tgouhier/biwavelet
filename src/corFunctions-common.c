/* 
 * Common functions for fast calculations of correlations
 *
 * General notes about handling missing data, zero MAD etc:
 * The idea is that bicor should switch to cor whenever it is feasible, it helps, and it is requested:
 * (1) if median is NA, the mean would be NA as well, so there's no point in switching to Pearson
 * (2) In the results, columns and rows corresponding to input with NA means/medians are NA'd out.
 * (3) The convention is that if zeroMAD is set to non-zero, it is the index of the column in which MAD is
 *     zero (plus one for C indexing)
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/libextern.h>
#define LDOUBLE 	long double

#include "pivot.h"
#include "corFunctions-common.h"

#define RefUX	0.5

/**
 * Here I first put all NAs to the end, then call the pivot function to find
 * the appropriate quantile of the remaining (finite) entries.
 * 
 * q is the quantile: 1/2 will give exactly the median above.
 */
double quantile(double * x, size_t n, double q, int copy, int * err)
{
  double * xx;
  double res;

  if (copy)
  {
    if ( (xx=malloc(n * sizeof(double)))==NULL ) 
    {
      Rprintf("Memory allocation error in quantile(). Could not allocate %d kB.\n", 
              (int) (n * sizeof(double) / 1024 + 1));
      *err = 1;
      return NA_REAL;
    }
    memcpy((void *)xx, (void *)x, n * sizeof(double));
  } else xx = x;

    
  *err = 0;
  // Put all NA's at the end.
  size_t bound = n;
  for (size_t i=n; i>0; ) 
  {
    i--;
    if (ISNAN(xx[i]))
    {
       bound--;
       xx[i] = xx[bound];
       xx[bound] = NA_REAL;
    }
  }

  // Rprintf("Quantile: q: %f, n: %d, bound: %d\n", q, n, bound);
  // Any non-NA's left?

  if (bound==0)
    res = NA_REAL;
  else
  // yes, return the appropriate pivot. 
    res = pivot(xx, bound, ( 1.0 * (bound-1))*q);

  if (copy) free(xx);

  return res;
}