// Function checks whether a given number is prime: returns 1 if it is prime, 0 otherwise

int IsPrime(ZZ n){
    ZZ s;
    ZZ t;
    t = 4;

    for(s = 2; t <= n; ++s){
        if(n % s == 0){ // n has factors other than n and 1, hence is composite.
            return(0);
        }
        else{
            t = t + 2 * s - 1;
        }
    }

    return(1); // n is prime no divisor other than n and 1 found
}