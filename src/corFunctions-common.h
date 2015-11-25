/* 
 * Common functions for fast calculations of correlations
 */
#ifndef __corFunctions_common_h__
#define __corFunctions_common_h__
#include "pivot.h"

enum { noWarning, warnZeroMAD };
double quantile(double * x, const size_t n, const double q, int * err);
#endif
