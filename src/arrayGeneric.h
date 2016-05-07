// The name of this file ends in .h so R CMD install does not attempt to compile
// it on its own.

class CLASS_NAME {

  protected:

    TYPE*   data_;
    size_t  size_;
    int 	  allocated;
    vector <size_t> dims;
    string	name_;

  public:

  #ifdef CheckDimensions

    TYPE value(const size_t i) {
      // fail fast
      if(i >= dims[1]) {
        stop("Index out of range in variable %s", name_);
      }

      return data_[i];
    }

    TYPE value(const size_t i, const size_t j) {
      // fail fast
      if(dims.size() != 2) {
        stop("incorrect number of dimensions accessing variable %s", name_);
      }

      // fail fast
      if ( (i >= dims[0]) || (j >= dims[1]) ) {
        stop("Index out of range in variable %s", name_);
      }

      return data_[j * dims[0] + i];
    }

    TYPE linValue(const size_t i) {
      size_t max = 1;

      for (size_t di = 0; di < dims.size(); di++) {
        max *= dims[di];
      }

      // fail fast
      if (i >= max) {
        stop("Linear index out of range in variable %s", name_);
      }

      return data_[i];
    }

    void setValue(size_t i, TYPE r) {
      if (i < dims[0]) {
        data_[i] = r;
      } else {
        stop("Index out of range in variable %s", name_);
      }
    }

    void setValue(size_t i, size_t j, TYPE r) {
      if (dims.size() == 2) {
        if ( (i < dims[0]) && (j < dims[1]) ) {
          data_[j * dims[0] + i] = r;
        } else {
          stop("Index out of range in variable %s", name_);
        }
      } else {
        stop("incorrect number of dimensions accessing variable %s", name_);
      }
    }

    void setValue(size_t i, size_t j, size_t k, TYPE r) {
      if (dims.size() == 3) {
        if ((k < dims[2]) && (j < dims[1]) && (i < dims[0])) {
          data_[(k * dims[1] + j) * dims[0] + i] = r;
        } else {
          stop("Index out of range in variable", name_);
        }
      } else {
        stop("Incorrect number of dimensions accessing variable %s", name_);
      }
    }

    void linValue(size_t i, TYPE r) {
      size_t max = 1;
      for (size_t di = 0; di < dims.size(); di++) {
        max *= dims[di];
      }

      if (i < max) {
        data_[i] = r;
      } else {
        stop("Linear index out of range in variable %s", name_);
      }
    }

  #else

    TYPE linValue(size_t i)
       { return data_[i]; }
    TYPE value(size_t i)
       { return data_[i]; }
    TYPE value(size_t i, size_t j)
       { return (data_[j * dims[0] + i]); }
    TYPE value(size_t i, size_t j, size_t k)
       { return (data_[(k * dims[1] + j) * dims[0] + i]); }

    void linValue(size_t i, TYPE r)
       { data_[i] = r; }
    void value(size_t i, TYPE r)
       { data_[i]; }
    void value(size_t i, size_t j, TYPE r)
       { data_[j * dims[0] + i] = r; }
    void value(size_t i, size_t j, size_t k, TYPE r)
       { data_[(k * dims[1] + j) * dims[0] + i] = r; }


  #endif

    void name(string n) {
      name_ = n;
    }

    string name() {
      return name_;
    }

    void setDim(size_t length);
    void setDim(size_t nrow, size_t ncol);
    void setDim(vector <size_t> dims, size_t start = 0);

    size_t nDim() {
      return dims.size();
    }

    vector <size_t> dim() {
      return dims;
    }

    size_t size() {
      return size_;
    }

    // Points the data_ pointer to given data.
    // Will make sure the data will not be deallocated in the destructor.
    void wrap(TYPE* data, size_t len)	{

      if (allocated) {
        delete data_;
      }

      allocated = 0;
      data_ = data;
      size_ = len;
      setDim(len);
    }

    void wrap(TYPE* data, size_t nrow, size_t ncol) {
      wrap(data, nrow * ncol);
      setDim(nrow, ncol);
    }

    TYPE* data() {
      return data_;
    }

    void copy2vector(size_t start, size_t length, vector<int>& result);
    void copy2vector(size_t start, size_t length, vector<double>& result);
    void rowQuantile(double q, dArray& quantile);

    CLASS_NAME() {
      allocated = 0;
      data_ = (TYPE*) NULL;
      dims.clear();
    }

    ~CLASS_NAME() {
      if (allocated) {
        delete data_;
        allocated = 0;
      }
    }

};

void CLASS_NAME::setDim(size_t length) {
  if (length > size_) {
    stop("attempt to set linear dimension %d higher than size %d in variable %s",
         length, size(), name());
  } else {
    dims.clear();
    dims.push_back(length);
  }
}

void CLASS_NAME::setDim(size_t nrow, size_t ncol) {
  if (nrow*ncol > size()) {
    stop("attempt to set matrix dimensions %d, %d higher than size %d in variable %s",
         nrow, ncol, size(), name());
  } else {
    dims.clear();
    dims.push_back(nrow);
    dims.push_back(ncol);
  }
}

void CLASS_NAME::setDim(vector <size_t> dims, size_t start) {
  size_t len = 1;
  for (size_t i = start; i < dims.size(); i++) {
    len *= dims[i];
  }

  if (len > size()) {
    stop("setDim: not enough space to accomodate given dimensions.");
  }

  this->dims.clear();
  this->dims.reserve(dims.size() - start);

  for (size_t i = start; i < dims.size(); i++) {
    this->dims.push_back(dims[i]);
  }
}
