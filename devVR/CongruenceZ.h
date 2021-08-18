
#include <math.h>
#include <vector>
#include <algorithm>

#include "NTL/ZZ.h"
#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
NTL_CLIENT

long powMod(long a, long n, long b) {
    // calculates a^n in O(log n)

    long ans = 1;            // Initialise answer

    while (n > 0) {
        if (n % 2 == 1) {   // if (n is odd) then
            ans = (ans * a) % b;
        }

        a = (a * a) % b;    // a = a^2 (mod b)
        n /= 2;             // n = n/2
    }

    return ans;
}

ZZ getMaxCoeff(ZZX x) {
    ZZ x_i = conv<ZZ>("0");
    for (long i = 0; i = deg(x); ++i) {
        if (x_i > coeff(x,i)) {
            x_i = coeff(x,i);
        }
        else {
            continue;
        }
    }

    return x_i;
}

ZZ evaluate(ZZX f, ZZ x) {
    // evaluates f(x) for x = x

    long fDeg = deg(f);
    ZZ ans = ConstTerm(f);

    for (long i = 1; i <= fDeg; ++i) {
        ans += coeff(f,i)*x;
    }

    return ans;
}

ZZX polyMultiply(ZZX f, ZZX g) {
    // polynomial multiplication by Binary Segmentation

    ZZ maxTermFG = to_ZZ(std::max(deg(f)+1,deg(g)+1));
    ZZ maxCoeffF = getMaxCoeff(f);
    ZZ maxCoeffG = getMaxCoeff(g);

    ZZ rhs = maxTermFG * maxCoeffF * maxCoeffG;         // rhs = max(U,V) * max(f_j) * max(g_k)

    ZZ b = conv<ZZ>("1");
    ZZ lhs = 2^b - 1;                                   // lhs = 2^b-1
    while (lhs < rhs) {                                 // Choose b such that 2^b − 1 > max(U,V) * max(x_i) * max(y_k)
        b = b + 1;
        lhs = 2^b - 1;
    }

    ZZ X = lhs;
    ZZ F = evaluate(f,X);                               // evaluate f(x) for x = 2^b - 1
    ZZ G = evaluate(g,X);                               // evaluate g(x) for x = 2^b - 1

    ZZ m = F * G;                                       // Integer multiply

    long fgDeg = deg(f) + deg(g) + 1;                   // Degree of polynomial product
    ZZ s[fgDeg];                                        // Store coefficients and constant in an array
    // std::vector<ZZ> s;                                  // Store coefficients and constant in vector
    // s.resize(fgDeg);

    for (long i = 0; i < fgDeg; ++i) {                  // Reassemble coefficients into signal
        s[i] = (m/(lhs^i)) % lhs;                       // Extract next b bits: s_i = floor( m/(2b−1)^i ) mod 2^b − 1
    }                                                   // N.B. -- NTL "/" operator floors result by default

    ZZX polyProduct;                                    // Base-b digits of m are desired coefficients
    for (long i = 0; i < fgDeg; ++i) {
        SetCoeff(polyProduct,i,s[fgDeg-i]);
    }

    polyProduct.normalize();                            // remove leading zeros on coefficients

    return polyProduct;
}

int CongruenceZ(const ZZ& n, const ZZ& r, const ZZ& r2) {
    // congruence test of polynomials in large integer form

    return 0;

}