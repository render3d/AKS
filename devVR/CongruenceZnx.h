
#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
NTL_CLIENT

int CongruenceZnx(const ZZ& n, const ZZ& r, const ZZ& r2){
    // congruence test of polynomials in regular form

    ZZ_p::init(n);                      //mod n
    ZZ_pX b = ZZ_pX(to_long(r), 1) - 1; // b = x^r - 1;
    ZZ_pX e = ZZ_pX(1, 1);              // e = x
    ZZ_pX d = PowerMod(e, n, b);        // d = x^n mod b, n

    for(long a = 1; a <= to_long(r2 - 1); ++a){

        ZZ_pX c = ZZ_pX(1, 1) - a ;         // c = x - a;
        ZZ_pX f = PowerMod(c, n, b);        // f =(x - a)^n mod c, n which is the RHS
        ZZ_pX g = d - a;                    // g = x^n - a mod c, n.

        if(f == g){
            return(1); // n is prime
        }
        else{
            return(a); // n is not prime.
        }
    }
}