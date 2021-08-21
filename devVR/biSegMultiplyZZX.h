
#include <math.h>
#include <vector>
#include <algorithm>

#include "NTL/ZZ.h"
#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
NTL_CLIENT

ZZ powMod(ZZ a, ZZ n, const ZZ& b) {
    // calculates a^n (mod b) in O(log n)
    // Tested on longs but not ZZ (to be NTL compatible)

    ZZ ans = ZZ(1);         // Initialise answer

    while (n > 0) {
        if (n % 2 == 1) {   // if (n is odd) then
            ans = (ans * a) % b;
        }

        a = (a * a) % b;    // a = a^2 (mod b)
        n /= 2;             // n = n/2
    }

    return ans;
}

ZZ getMaxCoeff(const ZZX& x) {
    // returns the largest coefficient of the input polynomial

    ZZ x_i = ConstTerm(x);

    for (long i = 1; i <= deg(x); ++i) {
        if (x_i < coeff(x,i)) {
            x_i = coeff(x,i);
        }
        else {
            continue;
        }
    }

    return x_i;
}

ZZ evaluate(const ZZX& f, const ZZ& x) {
    // evaluates f(y) for y = x

    long fDeg = deg(f);
    ZZ ans = ConstTerm(f);

    for (long i = 1; i <= fDeg; ++i) {
        ans += coeff(f,i)*power(x,i);
    }

    return ans;
}

ZZX polyMultiply(const ZZX& f, const ZZX& g) {
    // polynomial multiplication by Binary Segmentation

    ZZ termsF = to_ZZ(deg(f)+1);
    std::cout << "\nDegree of f is = " << deg(f) <<"\n";
    std::cout << "No. of terms in f = " << termsF <<"\n";

    ZZ termsG = to_ZZ(deg(g)+1);
    std::cout << "Degree of g is = " << deg(g) <<"\n";
    std::cout << "No. of terms in g = " << termsG <<"\n";

    ZZ maxTermFG = to_ZZ(std::max(termsF,termsG));
    std::cout << "\nmax(U,V) = " << maxTermFG <<"\n";

    ZZ maxCoeffF = getMaxCoeff(f);
    std::cout << "max(f_j) = " << maxCoeffF <<"\n";

    ZZ maxCoeffG = getMaxCoeff(g);
    std::cout << "max(g_k) = " << maxCoeffG <<"\n";

    ZZ rhs = maxTermFG * maxCoeffF * maxCoeffG;         // rhs = max(U,V) * max(f_j) * max(g_k)
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

    ZZX polyProduct;                                    // Base-b digits of m are desired coefficients
    polyProduct.SetLength(fgTrm);                       // set the length of the underlying coefficient vector to number of Terms in polynomial product
    std::cout << "\nPolynomial product preallocated: " << polyProduct << "\n";

    for (long j = 0; j < fgTrm; ++j) {
        std::cout << "Coefficient of x^i term when i = " << j <<": " << s[j] << "\n";
        SetCoeff(polyProduct,j,s[j]);
        std::cout << "Polynomial product after i = " << j <<": " << polyProduct << "\n\n";
    }

    polyProduct.normalize();                            // remove leading zeros on coefficients

    return polyProduct;
}

ZZX polyMulMod(const ZZX& f, const ZZX& g, const ZZX& h) {
    // polynomial multiply (mod h) by Binary Segmentation and where h is also a polynomial

    ZZ termsF = to_ZZ(deg(f)+1);
    ZZ termsG = to_ZZ(deg(g)+1);

    ZZ maxTermFG = to_ZZ(std::max(termsF,termsG));
    std::cout << "\nmax(U,V) = " << maxTermFG <<"\n";

    ZZ maxCoeffF = getMaxCoeff(f);
    std::cout << "max(f_j) = " << maxCoeffF <<"\n";

    ZZ maxCoeffG = getMaxCoeff(g);
    std::cout << "max(g_k) = " << maxCoeffG <<"\n";

    ZZ rhs = maxTermFG * maxCoeffF * maxCoeffG;         // rhs = max(U,V) * max(f_j) * max(g_k)
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

    ZZX polyProduct;                                    // Base-b digits of m are desired coefficients
    polyProduct.SetLength(fgTrm);                       // set the length of the underlying coefficient vector to number of Terms in polynomial product
    std::cout << "\nPolynomial product preallocated: " << polyProduct << "\n";

    for (long j = 0; j < fgTrm; ++j) {
        std::cout << "Coefficient of x^i term when i = " << j <<": " << s[j] << "\n";
        SetCoeff(polyProduct,j,s[j]);
        std::cout << "Polynomial product after i = " << j <<": " << polyProduct << "\n\n";
    }

    polyProduct.normalize();                            // remove leading zeros on coefficients

    ZZX polyProdMod = polyProduct % h;

    return polyProdMod;
}

ZZX polyMulModN(const ZZX& f, const ZZX& g, const ZZ& n) {
    // polynomial multiply (mod h) by Binary Segmentation and where h is also a polynomial

    ZZ termsF = to_ZZ(deg(f)+1);
    ZZ termsG = to_ZZ(deg(g)+1);

    ZZ maxTermFG = to_ZZ(std::max(termsF,termsG));
    std::cout << "\nmax(U,V) = " << maxTermFG <<"\n";

    ZZ maxCoeffF = getMaxCoeff(f);
    std::cout << "max(f_j) = " << maxCoeffF <<"\n";

    ZZ maxCoeffG = getMaxCoeff(g);
    std::cout << "max(g_k) = " << maxCoeffG <<"\n";

    ZZ rhs = maxTermFG * maxCoeffF * maxCoeffG;         // rhs = max(U,V) * max(f_j) * max(g_k)
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

    ZZX polyProdModN;                                    // Base-b digits of m are desired coefficients
    polyProdModN.SetLength(fgTrm);                       // set the length of the underlying coefficient vector to number of Terms in polynomial product
    std::cout << "\nPolynomial product preallocated: " << polyProdModN << "\n";

    for (long j = 0; j < fgTrm; ++j) {
        std::cout << "Coefficient of x^i term when i = " << j <<": " << s[j] << " mod n = " << (s[j] % n) << "\n";
        SetCoeff(polyProdModN,j,(s[j] % n));            // s[j] (mod n)
        std::cout << "Polynomial product after i = " << j <<": " << polyProdModN << "\n\n";
    }

    polyProdModN.normalize();                            // remove leading zeros on coefficients

    return polyProdModN;
}

ZZX polyPowMod(ZZX a, ZZ n, ZZX b) {
    // calculates a^n (mod b) in O(log n)
    // Tested on longs but types need modifying to be NTL compatible

    ZZX ans;                            // Initialise answer
    SetCoeff(ans,0,1);

    while (n > 0) {
        if (n % 2 == 1) {               // if (n is odd) then
            std::cout << "\n" << ans << "*" << a << " % " << b << " =\n";
            ans = polyMultiply(ans,a) % b;
            // ans = polyMulMod(ans,a,b);
            // ans %= b;
            std::cout << "\n" << ans << "\n";
        }
        std::cout << "\n" << a << "*" << a << " % " << b << " =\n";
        // a = polyMultiply(a,a) % b;    // a = a^2 (mod b)
        a = polyMultiply(a,a);
        // a %= b;
        std::cout << "\n" << a << "\n";
        n /= 2;                         // n = n/2
    }

    return ans;
}

ZZX polyPowModN(ZZX a, ZZ e, ZZ n) {

    ZZX ans;

    long aDeg = deg(a);
    for (long i = 0; i <= aDeg; ++i) {
        SetCoeff(ans,i,(coeff(a,i) % n));
    }

    return ans;
}

ZZX PowModNpX(ZZX a, ZZ e, ZZX b, ZZ n) {
    // calculates a^n (mod b, c) in O(log n)
    // Tested on longs but types need modifying to be NTL compatible

    ZZX ans;                            // Initialise answer
    SetCoeff(ans,0,1);

    while (e > 0) {
        std::cout << "\n" << "e = " << e << "\n";
        if (e % 2 == 1) {               // if (e is odd) then
            std::cout << "\n" << ans << "*" << a << " % " << b << " =\n";
            // ans = polyMultiply(ans,a) % b;
            ans = polyMulModN(ans,a,n);
            std::cout << "\n ans (mod n) = " << ans << "\n";
            ans = polyPowModN((ans % b),ZZ(1),n);
            std::cout << "\n ans (mod b, n)= " << ans << "\n";
        }
        std::cout << "\n" << a << "*" << a << " % " << b << " =\n";
        // a = polyMultiply(a,a) % b;    // a = a^2 (mod b)
        a = polyMulModN(a,a,n) % b;
        // a %= b;
        std::cout << "\n a (mod b, n)= " << a << "\n";
        e /= 2;                         // e = e/2
    }

    return ans;
}
