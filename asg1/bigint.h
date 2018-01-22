#ifndef __BIGINT_H__
#define __BIGINT_H__

#include <exception>
#include <iostream>
#include <utility>
using namespace std;

#include "debug.h"
using digitFunc = unsigned char;
using bigintFunc = vector<digitFunc>;

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
    friend void multiply_by_2 (unumber&);
    friend void divide_by_2 (unumber&);
    friend bigintFunc bigAdd (const bigintFunc&, const bigintFunc&);
    friend bigintFunc bigSub (const bigintFunc&, const bigintFunc&);
    friend bigintFunc bigMult (const bigintFunc&, const bigintFunc&);
    friend bool bigDiff (const bigintFunc&, const bigintFunc&);
    friend bigintFunc pprod(const bigintFunc&, size_t);
    friend bigintFunc pquot(const bigintFunc&, size_t);
    friend bigintFunc prem(const bigintFunc&, size_t);
    friend digitFunc testDigit(const bigintFunc&, const bigintFunc&, size_t, size_t);
    friend bool minFinder(const bigintFunc&, const bigintFunc&, size_t, size_t);
    friend bigintFunc diff(const bigintFunc&, const bigintFunc&, size_t, size_t);
    friend quot_rem longdiv(const bigintFunc& x, const bigintFunc& y, size_t n, size_t m);
    friend quot_rem div(const bigintFunc& x, const bigintFunc& y);
  
  
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
