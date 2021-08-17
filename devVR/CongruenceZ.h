
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

int CongruenceZ(const ZZ& n, const ZZ& r, const ZZ& r2) {
    // congruence test of polynomials in large integer form

    return 0;

}