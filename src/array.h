// The name of this file ends in h so R CMD install doesn't compile it twice. Not very clean but works
// for now.

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <R.h>
#include "Exception.h"

extern "C" {
  #include "quantile.h"
}

#ifndef __array_cc__
#define __array_cc__

// #define ISNAN(x) false
using namespace std;

// define a class that can conveniently hold a big vector, matrix etc.
#define NoDim	-1
#define CheckDimensions

string NumberToString(int n) {
  char s[100];
  string ss;
  sprintf(s, "%d", n);
  ss = s;
  return ss;
} 

class indArray {
  protected:

    size_t 	* data_;
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
        throw(Exception(string(
            "indArray::value: index out of range in variable") + name()));
      }

      const size_t j = (i % (8 * sizeof(size_t)));
      return (data_[ii] & mask[j]) != 0;
    }

    void value(const size_t i, const bool v) {
      const size_t ii = (i / (8 * sizeof(size_t)));
      if (ii >= size_) {
        throw(Exception(string(
            "indArray::value: index out of range in variable") + name()));
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

    void init(const size_t size);
    void init(const size_t size, const bool value);
    size_t size() { return size_ * 8 * sizeof(size_t); }
    size_t * data() { return data_; }

    void show() {
     cout << "data_:";
     for (size_t i = 0; i < size_; i++) {
       cout << data_[i] << ", ";
     }
     cout << endl;
   }

    indArray() {
      allocated = 0;
      data_ = (size_t *) NULL;
    }

    indArray(size_t size) { init(size); }

    indArray(size_t size, bool value) { init(size, value); }

    // indArray(indArray arr);      // This constructor will copy the data from arr into *this

    ~indArray() { if (allocated) { delete data_; allocated = 0; } }

};

void indArray::init(size_t size) {
  size_t val = 1;
  for (size_t i=0; i<8*sizeof(size_t); i++) {
    mask[i] = val; 
    invMask[i] = ~mask[i];
    if (i < 8*sizeof(size_t)-1) val*=2;
  }
  val32[0] = 0;
  val32[1] = ~0;
  size_ = (size_t)(ceil((size * 1.0) / (8 * sizeof(size_t))));
  data_ = new size_t[size_];
  allocated = 1;
}

void indArray::init(size_t size, bool value) {
  init(size);
  for (size_t i = 0; i < size_; i++) {
    data_[i] = val32[(size_t) value];
  }
}

#define INT_CLASS iArray

class dArray;

#define TYPE int
#define CLASS_NAME iArray
#include "arrayGeneric.h"
#undef TYPE
#undef CLASS_NAME

#define TYPE double
#define CLASS_NAME dArray
#include "arrayGeneric.h"
#undef TYPE
#undef CLASS_NAME

/**
 * Find maximum in a vector of doubles.
 */
double max(const vector <double> v) {
  if (v.size() == 0) {
    throw(Exception(string("attempt to calculate max of an empty vector.")));
  }

  double max = v[0];
  for (size_t i = 1; i < v.size(); i++) {
    if (!ISNAN(v[i]) && (v[i] > max)) {
      max = v[i];
    }
  }

  return max;
}

/**
 * Find maximum in a vector of ints.
 */
int max(const vector <int> v) {
  if (v.size()==0) {
    throw(Exception(string("attempt to calculate max of an empty vector.")));
  }

  int max = v[0];
  for (size_t i=1; i<v.size(); i++) {
    if (v[i] > max) {
      max = v[i];
    }
  }

  return max;
}

/**
 * Find minimum in a vector of doubles.
 */
double min(const vector <double> v) {
  if (v.size() == 0) {
    throw(Exception(string("attempt to calculate min of an empty vector.")));
  }
    
  double min = v[0];
  for (size_t i = 1; i < v.size(); i++) {
    if (!ISNAN(v[i]) && (v[i] < min)) {
      min = v[i];
    }
  }
    
  return min;
}

/**
 * Find minimum in a vector of ints.
 */
int min(const vector <int> v) {
  if (v.size() == 0) {
    throw(Exception(string("attempt to calculate min of an empty vector.")));
  }
  int min = v[0];
  for (size_t i = 1; i < v.size(); i++) {
    if (v[i] < min) min = v[i];
  }
  return min;
}

/**
 * Columnwise quantile of array.
 */
void dArray::colQuantile(const double q, dArray & quant) {

  if (dim().size() == 0) {
    throw(Exception(string(
       "Attempt to calculate columnwise quantile of array that has no dimensions set.")));
  }
  
  if (dim().size() == 1) {
    quant.setDim(1);
  } else {
    quant.setDim(dim(), 1);
  }
    

  const size_t colLen = dim()[0], totLen = length();

  if (colLen == 0) {
    throw(Exception(string(
        "colQuantile: Column length is zero in variable") + name()));
  }

  vector <double> column;
  column.reserve(colLen);

  int err;
  double val;
  
  for (size_t i = 0, col = 0; i < totLen; i += colLen, col++) {
    copy2vector(i, colLen, column);
    val = quantile(&(column[0]), colLen, q,  & err);
    quant.linValue(col, val);
  }
}

/**
 * Row-wise quantile of an array.
 */
void dArray::rowQuantile(const double q, dArray &quant) {
  if (dim().size() == 0) {
    throw(Exception(string(
        "Attempt to calculate row-wise quantile of array that has no dimensions set.")));
  }
  
  if (dim().size() == 1) {
    quant.setDim(1);
  } else {
    if (dim().size()>2) {
      throw(Exception(string(
       "Row-wise quantiles are only defined for 2-dimensional arrays.")));
    }
    
    vector <size_t> dim1 = dim();
    dim1.pop_back();
    quant.setDim(dim1);
  }

  const size_t rowLen = dim()[1];
  const size_t nrow = dim()[0];

  if (rowLen == 0) {
    throw(Exception(string(
        "rowQuantile: Row length is zero in variable") + name()));
  }

  vector <double> rowData;
  rowData.reserve(rowLen);

  int err;
  double val;
  for (size_t row = 0; row < nrow; row++) {
    rowData.clear();
    for (size_t col = 0; col < rowLen; col++) {
      rowData.push_back(value(row, col));
    }
    val = quantile(&(rowData[0]), rowLen, q, & err);
    quant.linValue(row, val);
  }
}

void iArray::colQuantile(const double q, dArray & quant) {
  if (dim().size() == 0) {
    throw(Exception(string(
        "Attempt to calculate columnwise quantile of array that has no dimensions set.")));
  }
  
  if (dim().size() == 1) {
    quant.setDim(1);
  } else {
    quant.setDim(dim(), 1);
  }

  const size_t colLen = dim()[0];
  const size_t totLen = length();

  if (colLen == 0) {
    throw(Exception(string(
        "colQuantile: Column length is zero in variable") + name()));
  }

  vector <double> column;
  column.reserve(colLen);

  int err;
  double val;
  for (size_t i = 0, col = 0; i < totLen; i += colLen, col++) {
    copy2vector(i, colLen, column);
    val = quantile(&(column[0]), colLen, q, & err);
    quant.linValue(col, val);
  }
}

#endif
