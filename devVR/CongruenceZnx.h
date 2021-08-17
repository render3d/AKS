
#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
// #include <NTL/BasicThreadPool.h>
NTL_CLIENT

unsigned int ncores = std::thread::hardware_concurrency(); // machine cores - may return 0 when not able to detect
// const auto SetNumThreads(ncores); // number of threads - should correspond to the number of available cores on your machine

long CongruenceZnx(const ZZ& n, const ZZ& r, const ZZ& r2, const long& a){
    // congruence test of polynomials in regular form

    ZZ_p::init(n);                      // initialise mod n

    ZZ_pX b = ZZ_pX(to_long(r), 1) - 1; // b = x^r - 1 (mod n);
    ZZ_pX e = ZZ_pX(1, 1);              // e = x (mod n)
    ZZ_pX d = PowerMod(e, n, b);        // d = x^n (mod b, n)

    // NTL_EXEC_RANGE(a,first,last)

        // for(long j = first; j <= last; ++j){
        for(long j = 1; j <= a; ++j){

            ZZ_pX c = ZZ_pX(1, 1) - j ;         // c = x - a (mod n);
            ZZ_pX f = PowerMod(c, n, b);        // f = (x - a)^n (mod b, n) - LHS
            ZZ_pX g = d - j;                    // g = x^n - a (mod b, n) - RHS

            if(f != g){
                return(to_long(j)); // n is not prime
            }
        }

    // NTL_EXEC_RANGE_END

    return(0); // n is prime
}