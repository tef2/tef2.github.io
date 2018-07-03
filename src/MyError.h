#ifndef MYERROR_H
#define MYERROR_H

#include <stdexcept>
#include <string>

using namespace std;

class MyError : public logic_error {
public:
  MyError(const string& msg = "") : logic_error(msg) {}
};

#endif
