#include <cstdlib>
#include <exception>
#include <limits>
#include <stack>
#include <stdexcept>
#include <algorithm>
using namespace std;

#include "bigint.h"
#include "debug.h"
#include "util.h"

#define LINE_LIMIT 69


bigint::bigint (const bigvalue_t& that): big_value (that) {
    DEBUGF ('~', this << " -> " << big_value << " (copy ctor)")
}

bigint::bigint (long that): long_value (that) {
    init(to_string(that));
}

bigint::bigint (const string& that) {
    init(that);
}



void bigint::init (const string& that) {
    negative = false;
    auto itor = that.cbegin();
    if (itor != that.cend()) {
        if (*itor == '_' || *itor == '-') {
            negative = true;
            ++itor;
        }
    }
    int newval = 0;
    for (; itor != that.end() ;) {
        newval = newval * 10 + *itor++ - '0';
    }
    auto ritor  = that.crbegin();
    while (ritor != that.crend() && *ritor != '_' ) {
        big_value.push_back(*ritor - '0');
        ritor++;
    }
    long_value = negative ? - newval : + newval;
}


bigvalue_t do_bigadd (const bigvalue_t& left, const bigvalue_t& right) {
    bigvalue_t sum;
    
    size_t min_size = min(left.size(), right.size());
    size_t i = 0;
    digit_t carry(0);
    digit_t digit_sum(0);
    while (i < min_size) {
        digit_sum = left.at(i) + right.at(i) + carry;
        if (digit_sum <= 9) {
            carry = 0;
        } else {
            carry = 1;
            digit_sum -= 10;
        }
        sum.push_back(digit_sum);
        i++;
    }
    for (; i < left.size(); i++) {
        digit_sum = left.at(i) + carry;
        if (digit_sum <= 9) {   
            carry = 0;
        } else {
            carry = 1;
            digit_sum -= 10;
        }
        sum.push_back(digit_sum);
    }
    for (; i < right.size(); i++) {
        digit_sum = right.at(i) + carry;
        if (digit_sum <= 9) {         
            carry = 0;
        } else {
              carry = 1;
            digit_sum -= 10;
        }
        sum.push_back(digit_sum);
        
    }
    if (carry != 1)
        return sum;
    else
        sum.push_back(1);

    return sum;
}

bigvalue_t do_bigsub (const bigvalue_t& left, const bigvalue_t& right) {
    bigvalue_t diff;
    digit_t borrow(0);
    digit_t digit_diff(0);
    size_t i = 0;
    
    while (i < right.size()) {
        
        if (left.at(i) - borrow >= right.at(i)) {
            digit_diff = left.at(i) - right.at(i) - borrow;
            borrow = 0;
        } else {
            digit_diff = 10 + left.at(i) - right.at(i) - borrow; 
            borrow = 1;
        }
        diff.push_back(digit_diff);
        i++;
    }
    for (; i < left.size(); i++) {
        if (left.at(i) >= borrow) {        
            digit_diff = left.at(i) - borrow;
            borrow = 0;
        } else {
            digit_diff = 10 + left.at(i) - borrow; 
            borrow = 1;
        }
        diff.push_back(digit_diff);
    }
  
        
    while (diff.size() > 1 && diff.back() == 0)
        diff.pop_back();
 
       
    return diff;
}


bool do_bigless (const bigvalue_t& left, const bigvalue_t& right) {
    if (left.size() > right.size())
        return false;
    else 
        return true;


    auto rit = right.crbegin();
    auto lit = left.crbegin();
    
    while (lit != left.crend()) {
        if (*lit > *rit)
            return false;
        else 
            return true;
        lit++; 
        rit++;
    }

    return false;
}

bigint operator+ (const bigint& left, const bigint& right) {
    bigint sum;
    if (left.negative != right.negative) {
        if (do_bigless(left.big_value, right.big_value)) {
            sum.big_value = do_bigsub(right.big_value, left.big_value);
            sum.negative = right.negative;
        } else { 
            sum.big_value = do_bigsub(left.big_value, right.big_value);
            sum.negative = left.negative;
        }
        return sum;
    } else {        
        sum.big_value = do_bigadd(left.big_value, right.big_value);
        sum.negative = left.negative;
        return sum;      
    }
}

bigint operator- (const bigint& left, const bigint& right) {
    bigint diff;
    if (left.negative == right.negative) {
        if (do_bigless(left.big_value, right.big_value)) {
            diff.big_value = do_bigsub(right.big_value, 
                    left.big_value);
            diff.negative = not right.negative;
        } else {
            diff.big_value = do_bigsub(left.big_value,
                    right.big_value);
            diff.negative = left.negative;
        }

        return diff;
    } else {
        diff.negative = left.negative;
        diff.big_value = do_bigadd(left.big_value, 
                right.big_value);
        return diff;
    }
}


bigint operator+ (const bigint& right) {
    bigint pos_bigint(right);
    pos_bigint.negative = false;
    return pos_bigint;
}

bigint operator- (const bigint& right) {
    bigint pos_bigint(right);
    pos_bigint.negative = true;
    return pos_bigint;
}

bigvalue_t do_bigmul (const bigvalue_t& left, const bigvalue_t& right) {
    bigvalue_t product(left.size() + right.size(), 0);
    digit_t c, d;
    for (size_t i = 0; i < left.size(); i++) {
        c = 0;
        for (size_t j = 0; j < right.size(); j++) {
            d = product.at(i+j) + (left.at(i) * right.at(j)) + c;
            product.at(i+j) = d % 10;
            c = d / 10;
        }
        product.at(i + right.size()) = c;
    }
    while (product.size() > 1 && product.back() == 0)
        product.pop_back();
    return product;
}

bigint operator* (const bigint& left, const bigint& right) {
    bigint product;
    product.negative = (left.negative != right.negative);
    product.big_value = do_bigmul(left.big_value, 
            right.big_value);
    return product;
}



bigvalue_t partial_prod(const bigvalue_t& x, size_t k) {
    int temp, carry;
    size_t size = x.size();
    bigvalue_t product(x.size() + 1, 0);
    carry = 0;
    for (size_t i = 0; i < size; i++) {
        temp = x.at(i) * k + carry;
        product.at(i) = temp % 10;
        carry = temp / 10;
    }
    product.at(size) = carry;
    while (product.size() > 1 && product.back() == 0)
        product.pop_back();
    DEBUGF ('/', "partial_prod(" << x << ", "
                 << (int) k << ") = " << product)
    return product;
}

bigvalue_t partial_quot(const bigvalue_t& x, size_t k) {
    int temp, carry;
    
    bigvalue_t quotient(x.size(), 0);
    carry = 0;
    size_t size = x.size();
    size_t i = size - 1;
    while (i < size) {
        temp = x.at(i) + 10 * carry;
        quotient.at(i) = temp / k;
        carry = temp % k;
        i--;
    }
    while (quotient.size() > 1 && quotient.back() == 0) {
        quotient.pop_back();
    }
    return quotient;
}

bigvalue_t partial_rem(const bigvalue_t& x, size_t k) {
    int carry = 0;
    size_t size = x.size();
    size_t i = size - 1;
    while (i < size) {
        carry = (x.at(i) + 10 * carry) % k;
        i--;
    }
    return bigvalue_t(1, carry);
}


digit_t trialdigit(const bigvalue_t& r, const bigvalue_t& d, size_t k, size_t m) {
    int d2, r3;
    size_t km = k + m;

    if ((r.size() <= km) && (r.size() <= km - 1) && (r.size() <= km - 2)) {
        r3 = 0;
    } else {
        if (r.size() > km) {
            r3 = (r.at(km)*10 + r.at(km - 1))*10 + r.at(km - 2);
        } else if (r.size() > km - 1) {
            r3 = r.at(km - 1)*10 + r.at(km - 2);
        } else if (r.size() > km - 2) {
            r3 = r.at(km - 2);
        }
    }

    if (d.size() > m - 1) {
        d2 = d.at(m - 1)*10 + d.at(m - 2);
    } else if (d.size() > m - 2) {
        d2 = d.at(m - 2);
    } else {
        d2 = 0;
    }

    return min(r3 / d2, 9);
}



bool smaller(const bigvalue_t& r, const bigvalue_t& dq, size_t k, size_t m) {
    int i, j;
    bigvalue_t r_copy(r);
    for (;r_copy.size() <= m + k;) {
        r_copy.push_back(0);
    }

    bigvalue_t dq_copy(dq);
    for (;dq_copy.size() <= m;) {
        dq_copy.push_back(0);
    }

    j = 0;
    i = m;
    
    while (i != j) {
        if (r_copy.at(i + k) == dq_copy.at(i)) {
            i = i - 1;           
        } else {
           j = i;
        } 
    }
    return r_copy.at(i + k) < dq_copy.at(i);
}

bigvalue_t difference(const bigvalue_t& r, const bigvalue_t& dq, size_t k, size_t m) {
    bigvalue_t dq_shifted;
    size_t i = 0;
    auto it = dq.cbegin();
    while (i < k) {
        dq_shifted.push_back(0);
        i++;
    }
    while (it != dq.cend()) {
        dq_shifted.push_back(*it);
        it++;
    }
    return do_bigsub (r, dq_shifted); 
}

bigint::quot_rem longdiv(const bigvalue_t& x, const bigvalue_t& y,size_t n, size_t m) {
    bigvalue_t d, dq, q(n, 0), r;
    int f, qt;
    int k;

    f = 10 / (y.at(m - 1) + 1);
    r = partial_prod(x, f);
    d = partial_prod(y, f);
    for (k = n - m; k >= 0; k--) {
        qt = trialdigit(r, d, k, m);
        dq = partial_prod(d, qt);
        if (smaller(r, dq, k, m)) {
            qt = qt - 1;
            dq = partial_prod(d, qt);
        }

        q.at(k) = qt;
        r = difference(r, dq, k, m);
    }

    while (r.size() > 1 && r.back() == 0)
        r.pop_back();
    while (q.size() > 1 && q.back() == 0)
        q.pop_back();
   
    return make_pair(bigint(q), bigint(partial_quot(r,f)));
}

bigint::quot_rem divide(const bigvalue_t& x, const bigvalue_t& y) {
    int n, y1; 
    int m = y.size();
    if (m != 1) {
        n = x.size();
        if (m > n)
            return make_pair(bigint(0), x);
        return longdiv(x, y, n, m);
    } else {
        y1 = y.at(m - 1);
        return make_pair(bigint(partial_quot(x, y1)), bigint(partial_rem(x, y1)));
    }
}

bigint operator/ (const bigint& left, const bigint& right) {
    if (right != 0) {
        bigint result = divide (left.big_value, right.big_value).first;
        result.negative = left.negative ^ right.negative;
        return result;
    } else {
        throw ydc_exn ("ydc: divid by zero");
    }  
}

bigint operator% (const bigint& left, const bigint& right) {
    if (left == NULL || right == NULL) {
        result.negative  = 0;
    } else { 
        bigint result = divide (left.big_value, right.big_value).second;
        result.negative = left.negative ^ right.negative;
    }
  
    return result;
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
    bool retval = false;
    if (left.negative) {
        if (right.negative) 
            retval = do_bigless(right.big_value, left.big_value);
        else 
            retval = true;
    } else {
        if (!right.negative)
            retval = do_bigless(left.big_value, right.big_value);           
        else 
            retval = false;
    }
    return retval;
}

ostream& operator<< (ostream& out, const bigint& that) {
    if (that.negative) {
         out << '-';
    } 
    out << that.big_value;
    return out;
}

ostream& operator<< (ostream& out, const bigvalue_t& that) {
    size_t i = 0;
    auto rit = that.crbegin();

    while (rit != that.crend()) {
        out << (int) *rit;
        i++;
        if (i == LINE_LIMIT) {
            out << "\\" << endl;
            i = 0;
        }
        ++rit;
    }
    return out;
}

long bigint::to_long() const {
    if (*this > bigint (numeric_limits<long>::min())) {
        if (*this <= bigint (numeric_limits<long>::max())) {
            return long_value;
        }
    } else {
        throw range_error ("bigint__to_long: out of range");
    }
    
    
}


bigint pow (const bigint& base, const bigint& exponent) {
    if (base != 0) {
        bigint base_copy = base;
        long expt = exponent.to_long();
        bigint result = 1;
        if (expt > 0) {
            for (; expt > 0 ; ) {
                if (expt & 1) { 
                    result = result * base_copy;
                    --expt;
                } else { 
                    base_copy = base_copy * base_copy;
                    expt /= 2;
                }
            }
        } else if (expt < 0) {
            base_copy = 1 / base_copy;
            expt = - expt;
        } 
        
        return result;
    } else {
        return 0;
    }
    
}
