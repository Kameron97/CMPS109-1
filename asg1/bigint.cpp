#include <cstdlib>
#include <exception>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <stack>

using namespace std;

#include "bigint.h"
#include "debug.h"
#include "util.h"



bigint::bigint (const bigintFunc& that): big_value (that) {
    DEBUGF ('~', this << " -> " << big_value << " (debugConstructor)")
}

bigint::bigint (long that): long_value (that) {
    init(to_string(that));
}

bigint::bigint (const string& that) {
    init(that);
}



void bigint::init (const string& that) {
    negative = false;
    auto i = that.cbegin();
    if (i != that.cend()) {
        if (*i == '_' || *i == '-') {
            negative = true;
            ++i;
        }
    }
    int value = 0;
    for (; i != that.end() ;) {
        value = value * 10 + *i++ - '0';
    }
    auto j  = that.crbegin();
    while (j != that.crend() && *j != '_' ) {
        big_value.push_back(*j - '0');
        j++;
    }
    long_value = negative ? - value : + value;
}


bigintFunc bigAdd (const bigintFunc& left, const bigintFunc& right) {
    bigintFunc sumAdd;
    size_t i = 0;
    digitFunc borrow(0);
    digitFunc digSum(0);
    while (i < min(left.size(), right.size())) {
        digSum = left.at(i) + right.at(i) + borrow;
        if (digSum <= 9) {
            borrow = 0;
        } else {
            borrow = 1;
            digSum -= 10;
        }
        sumAdd.push_back(digSum);
        i++;
    }
    for (; i < left.size(); i++) {
        digSum = left.at(i) + borrow;
        if (digSum <= 9) {   
            borrow = 0;
        } else {
            borrow = 1;
            digSum -= 10;
        }
        sumAdd.push_back(digSum);
    }
    for (; i < right.size(); i++) {
        digSum = right.at(i) + borrow;
        if (digSum <= 9) {         
            borrow = 0;
        } else {
            borrow = 1;
            digSum -= 10;
        }
        sumAdd.push_back(digSum);
        
    }
    if (borrow != 1)
        return sumAdd;
    else
        sumAdd.push_back(1);

    return sumAdd;
}

bigintFunc bigSub (const bigintFunc& left, const bigintFunc& right) {
    bigintFunc diffSub;
    digitFunc borrow(0);
    digitFunc digSub(0);
    size_t i = 0;
    
    while (i < right.size()) {
        
        if (left.at(i) - borrow >= right.at(i)) {
            digSub = left.at(i) - right.at(i) - borrow;
            borrow = 0;
        } else {
            digSub = 10 + left.at(i) - right.at(i) - borrow; 
            borrow = 1;
        }
        diffSub.push_back(digSub);
        i++;
    }
    for (; i < left.size(); i++) {
        if (left.at(i) >= borrow) {        
            digSub = left.at(i) - borrow;
            borrow = 0;
        } else {
            digSub = 10 + left.at(i) - borrow; 
            borrow = 1;
        }
        diffSub.push_back(digSub);
    }
  
        
    while (diffSub.size() > 1 && diffSub.back() == 0)
        diffSub.pop_back();
 
       
    return diffSub;
}


bool bigDiff (const bigintFunc& left, const bigintFunc& right) {
    if (left.size() > right.size()) {
        return false;
    } else {
        return true;
    }

    auto r = right.crbegin();
    auto l = left.crbegin();
    
    while (l != left.crend()) {
        if (*l > *r) {
            return false;
        } else { 
            return true;
        }
        l++; 
        r++;
    }

    return false;
}

bigint operator+ (const bigint& left, const bigint& right) {
    bigint bigSum;
    if (left.negative != right.negative) {
        if (bigDiff(left.big_value, right.big_value)) {
            bigSum.big_value = bigSub(right.big_value, left.big_value);
            bigSum.negative = right.negative;
        } else { 
            bigSum.big_value = bigSub(left.big_value, right.big_value);
            bigSum.negative = left.negative;
        }
        return bigSum;
    } else {        
        bigSum.big_value = bigAdd(left.big_value, right.big_value);
        bigSum.negative = left.negative;
        return bigSum;      
    }
}

bigint operator- (const bigint& left, const bigint& right) {
    bigint bigDiff;
    if (left.negative == right.negative) {
        if (bigDiff(left.big_value, right.big_value)) {
            bigDiff.big_value = bigSub(right.big_value, 
                    left.big_value);
            bigDiff.negative = not right.negative;
        } else {
            bigDiff.big_value = bigSub(left.big_value,
                    right.big_value);
            bigDiff.negative = left.negative;
        }

        return bigDiff;
    } else {
        bigDiff.negative = left.negative;
        bigDiff.big_value = bigAdd(left.big_value, 
                right.big_value);
        return bigDiff;
    }
}


bigint operator+ (const bigint& right) {
    bigint posPlus(right);
    posPlus.negative = false;
    return posPlus;
}

bigint operator- (const bigint& right) {
    bigint posDiff(right);
    posDiff.negative = true;
    return posDiff;
}

bigintFunc bigMult (const bigintFunc& left, const bigintFunc& right) {
    bigintFunc prod(left.size() + right.size(), 0);
    digitFunc count, dInc;
    for (size_t i = 0; i < left.size(); i++) {
        count = 0;
        for (size_t j = 0; j < right.size(); j++) {
            dInc = prod.at(i+j) + (left.at(i) * right.at(j)) + count;
            prod.at(i+j) = dInc % 10;
            count = dInc / 10;
        }
        prod.at(i + right.size()) = count;
    }
    while (prod.size() > 1 && prod.back() == 0)
        prod.pop_back();
    return prod;
}

bigint operator* (const bigint& left, const bigint& right) {
    bigint prod;
    prod.negative = (left.negative != right.negative);
    prod.big_value = bigMult(left.big_value, 
            right.big_value);
    return prod;
}



bigintFunc pprod(const bigintFunc& x, size_t k) {
    int carryOver = 0;
    bigintFunc prod(x.size() + 1, 0);
    int tmp = 0;
    for (size_t i = 0; i < x.size(); i++) {
        tmp = x.at(i) * k + carryOver;
        prod.at(i) = tmp % 10;
        carryOver = tmp / 10;
    }
    prod.at(x.size()) = carryOver;
    while (prod.size() > 1 && prod.back() == 0)
        prod.pop_back();
    return prod;
}

bigintFunc pquot(const bigintFunc& x, size_t k) {
    int carryOver = 0; 
    bigintFunc quotient(x.size(), 0);
    size_t i = x.size() - 1;
    int tmp = 0;
    while (i < x.size()) {
        tmp = x.at(i) + 10 * carryOver;
        quotient.at(i) = tmp / k;
        carryOver = tmp % k;
        i--;
    }
    while (quotient.size() > 1 && quotient.back() == 0) {
        quotient.pop_back();
    }
    return quotient;
}

bigintFunc prem(const bigintFunc& x, size_t k) {
    int carryOver = 0;
    size_t i = x.size() - 1;
    while (i < x.size()) {
        carryOver = (x.at(i) + 10 * carryOver) % k;
        i--;
    }
    return bigintFunc(1, carryOver);
}


digitFunc testDigit(const bigintFunc& r, const bigintFunc& d, size_t k, size_t m) {
    int rInc = 0;
    size_t kInc = k + m;

    if ((r.size() <= kInc) && (r.size() <= kInc - 1) && (r.size() <= kInc - 2)) {
        rInc = 0;
    } else {
        if (r.size() > kInc) {
            rInc = (r.at(kInc)*10 + r.at(kInc - 1))*10 + r.at(kInc - 2);
        } else if (r.size() > kInc - 1) {
            rInc = r.at(kInc - 1)*10 + r.at(kInc - 2);
        } else if (r.size() > kInc - 2) {
            rInc = r.at(kInc - 2);
        }
    }

    int dInc = 0;
    if (d.size() > m - 1) {
        dInc = d.at(m - 1)*10 + d.at(m - 2);
    } else if (d.size() > m - 2) {
        dInc = d.at(m - 2);
    } else {
        dInc = 0;
    }

    return min(rInc / dInc, 9);
}



bool minFinder(const bigintFunc& r, const bigintFunc& dq, size_t k, size_t m) {
    bigintFunc copy(r);
    for (;copy.size() <= m + k;) {
        copy.push_back(0);
    }

    bigintFunc dq_copy(dq);
    for (;dq_copy.size() <= m;) {
        dq_copy.push_back(0);
    }

    int j = 0;
    int i = m;
    
    while (i != j) {
        if (copy.at(i + k) == dq_copy.at(i)) {
            i = i - 1;           
        } else {
           j = i;
        } 
    }
    return copy.at(i + k) < dq_copy.at(i);
}

bigintFunc diff(const bigintFunc& r, const bigintFunc& dq, size_t k, size_t m) {
    bigintFunc shift;
    size_t i = 0;
    auto dqBegin = dq.cbegin();
    while (i < k) {
        shift.push_back(0);
        i++;
    }
    while (dqBegin != dq.cend()) {
        shift.push_back(*dqBegin);
        dqBegin++;
    }
    return bigSub (r, shift); 
}

bigint::quot_rem longdiv(const bigintFunc& x, const bigintFunc& y,size_t n, size_t m) {
    bigintFunc dInc, dInc2, qtr(n, 0), rInc;
    int fInc, qInc;
    int kInc;

    fInc = 10 / (y.at(m - 1) + 1);
    rInc = pprod(x, fInc);
    dInc = pprod(y, fInc);
    for (kInc = n - m; kInc >= 0; kInc--) {
        qInc = testDigit(rInc, dInc, kInc, m);
        dInc2 = pprod(dInc, qInc);
        if (minFinder(rInc, dInc2, kInc, m)) {
            qInc = qInc - 1;
            dInc2 = pprod(dInc, qInc);
        }

        qtr.at(kInc) = qInc;
        rInc = diff(rInc, dInc2, kInc, m);
    }

    while (rInc.size() > 1 && rInc.back() == 0)
        rInc.pop_back();
    while (qtr.size() > 1 && qtr.back() == 0)
        qtr.pop_back();
   
    return make_pair(bigint(qtr), bigint(pquot(rInc,fInc)));
}

bigint::quot_rem divide(const bigintFunc& x, const bigintFunc& y) {
    if (y.size() != 1) {
        if (y.size() > x.size()) {
            return make_pair(bigint(0), x);
        } else {
            return longdiv(x, y, x.size(), y.size());
        }       
    } else {
        int yInc = y.at(y.size() - 1);
        return make_pair(bigint(pquot(x, yInc)), bigint(prem(x, yInc)));
    }
}

bigint operator/ (const bigint& left, const bigint& right) {
    if (right != 0) {
        bigint res = divide (left.big_value, right.big_value).first;
        res.negative = left.negative ^ right.negative;
        return res;
    } else {
        throw ydc_exn ("ydc error: Cannot divide by zero");
    }  
}

bigint operator% (const bigint& left, const bigint& right) {
    bigint res = divide (left.big_value, right.big_value).second;
    res.negative = left.negative ^ right.negative;
    return res;
}

bool operator== (const bigint& left, const bigint& right) {
    if (left.negative == right.negative || left.big_value.size() == right.big_value.size()) {
        return false;
    } else {
        return true;
    }
    size_t i = 0;
    while (i < left.big_value.size()) {
        if (left.big_value.at(i) != right.big_value.at(i)) {
            return false;
        }
        i++;
    }

    return true;
}

bool operator< (const bigint& left, const bigint& right) {
    bool isCorrectOperator = false;
    if (left.negative) {
        if (right.negative) 
            isCorrectOperator = bigDiff(right.big_value, left.big_value);
        else 
            isCorrectOperator = true;
    } else {
        if (!right.negative)
            isCorrectOperator = bigDiff(left.big_value, right.big_value);           
        else 
            isCorrectOperator = false;
    }
    return isCorrectOperator;
}

ostream& operator<< (ostream& out, const bigint& that) {
    if (that.negative) {
         out << '-';
    } 
    out << that.big_value;
    return out;
}

ostream& operator<< (ostream& out, const bigintFunc& that) {
    size_t i = 0;
    auto r = that.crbegin();

    while (r != that.crend()) {
        out << (int) *r;
        i++;
        if (i == 69) {
            out << "\\" << endl;
            i = 0;
        }
        ++r;
    }
    return out;
}

long bigint::to_long() const {
    if (*this > bigint (numeric_limits<long>::min())) {
        if (*this <= bigint (numeric_limits<long>::max())) {
            return long_value;
        }
    } else {
        throw range_error ("Out of range");
    }
    
    
}


bigint pow (const bigint& base, const bigint& exponent) {
    if (base != 0) {
        bigint copy = base;
        long power = exponent.to_long();
        bigint res = 1;
        if (power > 0) {
            for (; power > 0 ; ) {
                if (power & 1) { 
                    res = res * copy;
                    --power;
                } else { 
                    copy = copy * copy;
                    power /= 2;
                }
            }
        } else if (power < 0) {
            copy = 1 / copy;
            power = - power;
        }        
        return res;
    } else {
        return 0;
    }
    
}