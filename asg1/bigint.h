/*
Aravind Patnam
apatnam@ucsc.edu
ASG1 = DC Calc
January 22nd, 2018
*/




#ifndef __BIGINT_H__
#define __BIGINT_H__

#include <exception>
#include <iostream>
#include <utility>
using namespace std;

#include "debug.h"
using digit_t = unsigned char;
using bigintFunc = vector<digit_t>;

class bigint {
    friend ostream& operator<< (ostream&, const bigint&);
  
  
  private:
    void init (const string&);
    long long_value {};
    bool negative;
    bigintFunc big_value;
    using quot_rem = pair<bigint,bigint>;
    using unumber = unsigned long;
    friend quot_rem divide (const bigint&, const bigint&);
    friend bigintFunc bigAdd (const bigintFunc&, const bigintFunc&);
    friend bigintFunc bigSub (const bigintFunc&, const bigintFunc&);
    friend bigintFunc bigMult (const bigintFunc&, const bigintFunc&);
    friend bool isLess (const bigintFunc&, const bigintFunc&);
    friend bigintFunc pprod(const bigintFunc&, sizeVector);
    friend bigintFunc pquot(const bigintFunc&, sizeVector);
    friend bigintFunc prem(const bigintFunc&, sizeVector);
    friend digit_t isCorrectDigit(const bigintFunc&, const bigintFunc&, sizeVector, sizeVector);
    friend bool testMin(const bigintFunc&, const bigintFunc&, sizeVector, sizeVector);
    friend bigintFunc findDiff(const bigintFunc&, const bigintFunc&, sizeVector, sizeVector);
    friend quot_rem longdiv(const bigintFunc& x, const bigintFunc& y,sizeVector n, sizeVector m);
    friend quot_rem divide(const bigintFunc& x, const bigintFunc& y);
  
  
  public:

    bigint() = default;
    bigint (const bigint&) = default;
    bigint (bigint&&) = default;
    bigint& operator= (const bigint&) = default;
    bigint& operator= (bigint&&) = default;
    ~bigint() = default;

    bigint (const bigintFunc&);
    bigint (const long);
    bigint (const string&);

    friend bigint operator+ (const bigint&, const bigint&);
    friend bigint operator- (const bigint&, const bigint&);
    friend bigint operator+ (const bigint&);
    friend bigint operator- (const bigint&);
    long to_long() const;


    friend bigint operator* (const bigint&, const bigint&);
    friend bigint operator/ (const bigint&, const bigint&);
    friend bigint operator% (const bigint&, const bigint&);

 
    friend bool operator== (const bigint&, const bigint&);
    friend bool operator<  (const bigint&, const bigint&);
};



ostream& operator<< (ostream& out, const bigintFunc& that);

bigint pow (const bigint& base, const bigint& exponent);

inline bool operator!= (const bigint &left, const bigint &right) {
    return not (left == right);
}
inline bool operator>  (const bigint &left, const bigint &right) {
    return right < left;
}
inline bool operator<= (const bigint &left, const bigint &right) {
    return not (right < left);
}
inline bool operator>= (const bigint &left, const bigint &right) {
    return not (left < right);
}
#endif
