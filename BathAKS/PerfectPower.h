#include <NTL/RR.h>

// function to calculate if n = aˆb
// takes input ZZ n and returns 1 if n is a perfect power, 0 otherwise

int PerfectPower(ZZ n){
    long b = 2;
    RR k = to RR (log(n) / log(2));
    ZZ a;

    while(b <= to long(k)){
        // b cannot be bigger than k
        long c = to_long(ceil((log(n) / log(2)) / to_long(b)));
        a = pow(2, c); // assign guess value for a

        while(power(a, b) > n){
            double d = to_double(((b - 1) * a + n / power(a , ( b-1))) / b); // Apply Integer Newton's Method
            ZZ e = to_ZZ(floor(d));
            a = to_long (e); // adjust a
        }

        if(n == power(a, b)){
            // if n is a perfect power.
            cout << n << " is a perfect power, hence is not prime \n" ;
            cout << "n = aˆb \ n" ;
            cout << "b = " << b << "\n" ;
            cout << "a = " << a << "\n\n" ;
            return ( 1 ) ;
        }
        else{
            b = b + 1;
        }
    }

    if(n != power(a,b)){
        return(0); // n is not a perfect power .
    }

}