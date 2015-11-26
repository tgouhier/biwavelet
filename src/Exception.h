#ifndef __Exception_h__
#define __Exception_h__

#include <iostream>

using namespace std;

/**
 * Exception handling. Just adding a bit more information to the standard
 * exception.
 */
class Exception {

  protected:
    string	_what;

  public:
    virtual string what() const throw() { return _what; }
    Exception(string wht) throw() {
      _what = wht;
    }
    ~Exception() throw() {}
};

#endif
