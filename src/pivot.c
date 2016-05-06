/*
 * Compilation:
 *  gcc --std=c99 -fPIC -O3 -o pivot.so -shared pivot.c
 */
#include <stdlib.h>
#include <math.h>
#include "pivot.h"

double vMax(double *v, const size_t len) {
  double mx = v[0];
  for (size_t i = 1; i < len; i++) {
    if (v[i] > mx) {
      mx = v[i];
    }
  }
  return mx;
}

double vMin(double * v, const size_t len) {
  double mn = v[0];
  for (size_t i = 1; i < len; i++) {
    if (v[i] < mn) {
      mn = v[i];
    }
  }
  return mn;
}

double pivot(double * v, const size_t len, const double target) {
  if (len > 2) {

    // pick the pivot, say as the median of the first, middle and last
    const size_t i1 = 0;
    const size_t i2 = len - 1;
    const size_t i3 = (len - 1) / 2;
    size_t ip;

    if (v[i1] <= v[i2]) {
      if (v[i2] <= v[i3])
        ip = i2;
      else if (v[i3] >= v[i1])
         ip = i3;
      else
         ip = i1;
    } else {
      if (v[i1] <= v[i3])
        ip = i1;
      else if (v[i2] <= v[i3])
        ip = i3;
      else
        ip = i2;
    }

    // put ip at the end
    double vp = v[ip];
    v[ip] = v[len - 1];
    v[len - 1] = vp;

    // pivot everything else
    size_t bound = 0;
    for (size_t i=0; i<len; i++) {
      if (v[i] < vp) {
        const double x = v[bound];
        v[bound] = v[i];
        v[i] = x;
        bound++;
      }
    }

    v[len - 1] = v[bound];
    v[bound] = vp;

    // Did we find the target?
    const double crit = target - bound;

    if (fabs(crit) > 1.0) {
      if (crit < 0) {
        return pivot(v, bound, target);
      } else {
        return pivot(v + bound + 1, len - bound - 1, target - bound - 1);
      }
    }

    if (crit < 0) {
       const double v1 = vMax(v, bound);
       return (v1 * (-crit) + vp * (1 + crit));
    } // else
    const double v2 = vMin(v + bound + 1, len - bound - 1);
    return (vp * (1 - crit) + v2 * crit);

  } else if (len == 2) {

      const double v1 = vMin(v, 2);
      const double v2 = vMax(v, 2);

      if (target < 0) {
        return v1;
      } else if (target > 1) {
        return v2;
      } else {
        return target * v2 + (1-target) * v1;
      }
  } // else
  return v[0];
}
