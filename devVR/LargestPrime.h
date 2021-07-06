// calculates the largest prime factor of an input n

ZZ LargestPrime(ZZ n){
    ZZ f; f = 1;
    ZZ r; r = 2;
    ZZ x = n ;

    while(x != 1 & sqr(r) <= n){
        while(x % r == 0){
            x = x / r;
            f = r;
        }
        r = r + 1;
    }

    if(x == 1){
        return (f) ;
    }
    else{
        return(x);
    }
}