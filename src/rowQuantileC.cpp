#include "array.h"
#include <R.h>

extern "C" {
  void rowQuantileC(double *data,
                    const int *nrow, const int *ncol,
                    const double *q, double *res) {
    
    const int nr = *nrow;
    const int nc = *ncol;
    
    dArray d;
    d.wrap(data, nr, nc);
    
    if ((*q<0) || (*q>1)) {
      Rf_error("quantileC: given quantile is out of range 0 to 1.");
      return;
    }
    
    dArray quant;
    quant.wrap(res, nr);
    
    d.rowQuantile(*q, quant);
  }
} // extern "C"
