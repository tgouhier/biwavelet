#ifndef __array_cc__
#define __array_cc__

#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

extern "C" {
#include "quantile.h"
}

// define a class that can conveniently hold a big vector, matrix etc.
#define NoDim	-1
#define CheckDimensions

class indArray {
  protected:

    size_t* data_;
    size_t 	size_;
    int 	  allocated;
    string	name_;

    size_t	val32[2];
    size_t	mask[8 * sizeof(size_t)];
    size_t	invMask[8 * sizeof(size_t)];

  public:

#ifdef CheckDimensions

    bool value(const size_t i) {
      const size_t ii = (i / (8 * sizeof(size_t)));
      if (ii >= size_) {
        stop("indArray::value: index out of range in variable %s", name());
      }

      const size_t j = (i % (8 * sizeof(size_t)));
      return (data_[ii] & mask[j]) != 0;
    }

    void value(const size_t i, const bool v) {
      const size_t ii = (i / (8 * sizeof(size_t)));
      if (ii >= size_) {
        stop("indArray::value: index out of range in variable %s", name());
      }

      const size_t j = (i % (8 * sizeof(size_t)));

      if (v) {
        data_[ii] |= mask[j];
      } else {
        data_[ii] &= invMask[j];
      }
    }

#else

    bool value(size_t i) {
      const size_t ii = (i / (8 * sizeof(size_t)));
      const size_t j = (i % (8 * sizeof(size_t)));
      return data_[ii] & mask[j] != 0;
    }

    void value(size_t i, bool v) {
      const size_t ii = (i / (8 * sizeof(size_t)));
      const size_t j = (i % (8 * sizeof(size_t)));
      if (v) {
        data_[ii] |= mask[j];
      } else {
        data_[ii] &= invMask[j];
      }
    }

#endif

    void name(string n) {
      name_ = n;
    }

    string name() {
      return name_;
    }

    size_t size() {
      return size_ * 8 * sizeof(size_t);
    }

    size_t * data() {
      return data_;
    }

    void show() {
     cout << "data_:";
     for (size_t i = 0; i < size_; i++) {
       cout << data_[i] << ", ";
     }
     cout << endl;
   }

    indArray() {
      allocated = 0;
      data_ = (size_t*) NULL;
    }

    ~indArray() {
      if (allocated) {
        delete data_;
        allocated = 0;
      }
    }
};

class dArray;

#define TYPE double
#define CLASS_NAME dArray
#include "arrayGeneric.h"
#undef TYPE
#undef CLASS_NAME

/**
 * Row-wise quantile of an array.
 */
void dArray::rowQuantile(const double q, dArray& quant) {
  if (dim().size() == 0) {
    stop("Attempt to calculate row-wise quantile of array that has no dimensions set.");
  }

  if (dim().size() == 1) {
    quant.setDim(1);
  } else {
    if (dim().size() > 2) {
      stop("Row-wise quantiles are only defined for 2-dimensional arrays.");
    }

    vector<size_t> dim1 = dim();
    dim1.pop_back();
    quant.setDim(dim1);
  }

  const size_t rowLen = dim()[1];
  const size_t nrow = dim()[0];

  if (rowLen == 0) {
    stop("rowQuantile: Row length is zero in variable %s", name());
  }

  vector<double> rowData;
  rowData.reserve(rowLen);

  int err;
  double val;
  for (size_t row = 0; row < nrow; row++) {
    rowData.clear();
    for (size_t col = 0; col < rowLen; col++) {
      rowData.push_back(value(row, col));
    }
    val = quantile(&(rowData[0]), rowLen, q, &err);
    quant.linValue(row, val);
  }
}

#endif
