
#include <math.h>
#include <vector>
#include <algorithm>

#include "NTL/ZZ.h"
#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
NTL_CLIENT

ZZ_p getMaxCoeff(const ZZ_pX& x) {
    // returns the largest coefficient of the input polynomial

    ZZ_p x_i = ConstTerm(x);
    ZZ comp = rep(ConstTerm(x));

    for (long i = 1; i <= deg(x); ++i) {
        if (comp < rep(coeff(x,i))) {
            x_i = coeff(x,i);
            comp = rep(coeff(x,i));
        }
        else {
            continue;
        }
    }

    return x_i;
}

ZZ evaluate(const ZZ_pX& f, const ZZ& x) {
    // evaluates f(y) for y = x

    long fDeg = deg(f);
    ZZ ans = rep(ConstTerm(f));
    std::cout << "Sum = " << ans << "\n";

    for (long i = 1; i <= fDeg; ++i) {
        ans += rep(coeff(f,i))*power(x,i);
        std::cout << "Sum = " << ans << "\n";
    }

    return ans;
}

ZZ_pX polyModMul(const ZZ_pX& f, const ZZ_pX& g) {
    // polynomial multiplication by Binary Segmentation

    ZZ termsF = to_ZZ(deg(f)+1);
    std::cout << "\nDegree of f is = " << deg(f) <<"\n";
    std::cout << "No. of terms in f = " << termsF <<"\n";

    ZZ termsG = to_ZZ(deg(g)+1);
    std::cout << "Degree of g is = " << deg(g) <<"\n";
    std::cout << "No. of terms in g = " << termsG <<"\n";

    ZZ maxTermFG = to_ZZ(std::max(termsF,termsG));
    std::cout << "\nmax(U,V) = " << maxTermFG <<"\n";

    ZZ_p maxCoeffF = getMaxCoeff(f);
    std::cout << "max(f_j) = " << maxCoeffF <<"\n";

    ZZ_p maxCoeffG = getMaxCoeff(g);
    std::cout << "max(g_k) = " << maxCoeffG <<"\n";

    ZZ rhs = maxTermFG * rep(maxCoeffF) * rep(maxCoeffG);         // rhs = max(U,V) * max(f_j) * max(g_k)
    std::cout << "\nRHS = " << rhs <<"\n";

    long b = 1;
    ZZ lhs = (power2_ZZ(b)) - 1;                        // lhs = 2^b-1
    while (lhs <= rhs) {                                 // Choose b such that 2^b − 1 > max(U,V) * max(x_i) * max(y_k)
        b = b + 1;
        lhs = (power2_ZZ(b)) - 1;
    }
    std::cout << "LHS = 2^" << b << " - 1 = " << lhs <<"\n\n";

    ZZ X = lhs;
    ZZ F = evaluate(f,X);                               // evaluate f(x) for x = 2^b - 1
    std::cout << "F = f(X) = f(" << X << ") = " << F <<"\n";
    ZZ G = evaluate(g,X);                               // evaluate g(x) for x = 2^b - 1
    std::cout << "G = g(X) = g(" << X << ") = " << G <<"\n\n";

    ZZ m = F * G;                                       // Integer multiply
    std::cout << "m = F * G = " << F << " * " << G << " = " << m <<"\n\n";

    // long fgDeg = deg(f) + deg(g);                       // Degree of polynomial product
    long fgTrm = deg(f) + deg(g) + 1;                   // Terms in polynomial product
    std::cout << "Terms in f(x) * g(x) = " << fgTrm << "\n";

    ZZ s[fgTrm];                                        // Store coefficients and constant in an array
    // std::vector<ZZ> s;                                  // Store coefficients and constant in vector
    // s.resize(fgTrm);

    // for (long i = 0; i < fgTrm; ++i) {                  // Reassemble coefficients into signal
    //     s[i] = (m/(power(lhs,i))) % lhs;                // Extract next b bits: s_i = floor( m/(2b−1)^i ) mod 2^b − 1
    //     std::cout << "Next " << b << " bits for i = " << i << " are " << s[i] << "\n";
    // }                                                   // N.B. -- NTL "/" operator floors result by default

    // s[0] = m % lhs;
    // for (long i = 1; i < fgTrm; ++i) {
    //     m = m / lhs;
    //     s[i] = m % lhs;
    // }

    s[0] = m % lhs;                                     // Reassemble coefficients into signal
    std::cout << "\nConstant is equal to " << s[0] << "\n";
    if (s[0] > lhs/2) {                                 // Extract next b bits: s_i = floor( m/(2b−1)^i ) mod 2^b − 1
        std::cout << "s[0] > (2^b - 1)/2 " << s[0] << " > " << lhs/2 << "\n";
        s[0] = s[0] - lhs;
        std::cout << "Therefore s[0] = s[0] - (2^b - 1) = " << s[0] << "\n";
        m = m + lhs;
        std::cout << "And m = m + (2^b - 1) = " << m << "\n";
    }
    for (long i = 1; i < fgTrm; ++i) {
        m = m / lhs;                                    // N.B. -- NTL "/" operator floors result by default
        std::cout << "For i = " << i << ", m = m / (2^b - 1) = " << m << "\n";

        s[i] = m % lhs;
        std::cout << "And s[i] = m % (2^b - 1) = " << s[i] << "\n";

        if (s[i] > lhs/2) {
            s[i] = s[i] - lhs;
            m = m + lhs;
        }
    }

    ZZ_pX polyModProd;                                    // Base-b digits of m are desired coefficients
    polyModProd.SetLength(fgTrm);                       // set the length of the underlying coefficient vector to number of Terms in polynomial product
    std::cout << "\nPolynomial product preallocated: " << polyModProd << "\n";

    for (long j = 0; j < fgTrm; ++j) {
        std::cout << "Coefficient of x^i term when i = " << j <<": " << s[j] << "\n";
        SetCoeff(polyModProd,j,to_ZZ_p(s[j]));
        std::cout << "Polynomial product after i = " << j <<": " << polyModProd << "\n\n";
    }

    polyModProd.normalize();                            // remove leading zeros on coefficients

    return polyModProd;
}

ZZ_pX ZZpPowMod(ZZ_pX a, ZZ n, ZZ_pX b) {
    // calculates a^n (mod b) in O(log n)
    // Tested on longs but types need modifying to be NTL compatible

    ZZ_pX ans;                            // Initialise answer
    SetCoeff(ans,0,1);

    while (n > 0) {
        if (n % 2 == 1) {               // if (n is odd) then
            std::cout << "\n" << ans << "*" << a << " % " << b << " =\n";
            ans = polyModMul(ans,a) % b;
            // ans %= b;
            std::cout << "\n" << ans << "\n";
        }
        std::cout << "\n" << a << "*" << a << " % " << b << " =\n";
        a = polyModMul(a,a) % b;    // a = a^2 (mod b)
        // a %= b;
        std::cout << "\n" << a << "\n";
        n /= 2;                         // n = n/2

    }

    return ans;
}
