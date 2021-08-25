
#include "NTL/ZZ.h"
#include "NTL/ZZ_p.h"
NTL_CLIENT

// making use of Fermats Theorem
// aˆlambda(r)=1 mod r for gcd(a,r) = 1;

ZZ EulerPhi(ZZ r){
    ZZ a = to_ZZ(2);
    ZZ_p::init(r); // work mod r

    int q=0;

    // want to find an a coprime to r
    while(GCD(r,a)!=1){
        a=a+1;
    }

    // a <= r by construction
    ZZ b =to_ZZ(1); // b is the order of a mod r.
    ZZ g; // g is aˆb mod r
    int s=0; // act as an indicator

    while(s==0){
        ZZ g1 = power_ZZ(to_long(a),to_long(b));
        ZZ_p g2 = to_ZZ_p(g1);

        if (g2==1){ // if aˆb=1 then b is phi(r) and
            //so we are done.
            s=1; // s=1 => break out of while loop
        }
        else{
            s=0; // else increase value of b.
            b=b+1;
        }
    }

    return(b); // b is phi(r), return this value to main.
}