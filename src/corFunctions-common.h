/* 
 * Common functions for fast calculations of correlations
 */
#ifndef __corFunctions_common_h__
#define __corFunctions_common_h__
#define LDOUBLE 	long double
#include "pivot.h"

enum { noWarning, warnZeroMAD };
double quantile(double * x, size_t n, double q, int copy, int * err);
#endif
