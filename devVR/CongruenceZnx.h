
#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
NTL_CLIENT

int CongruenceZnx(const ZZ& n, const ZZ& r, const ZZ& r2){
    // congruence test of polynomials in regular form

    for(long a = 1; a <= to_long(r2 - 1); ++a){
        ZZ_p::init(n); //mod n
        ZZ_pX b = ZZ_pX(to_long(r), 1) - 1; // b = x^r - 1;
        ZZ_pX c = ZZ_pX(1, 1) - a ; // c = x - a;
        ZZ_pX f = PowerMod(c, n, b); // f =(x - a)^n mod c, n which is the RHS
        ZZ_pX e = ZZ_pX(1, 1);
        ZZ_pX g = PowerMod(e, n, b); // x^n mod b, n
        g = g - a ; // g1 = x^n - a mod c, n.

        if(f == g){
            return(1); // n is prime
        }
        else{
            return(a); // n is not prime.
        }
    }
}