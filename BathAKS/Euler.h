// Euler's phi function

ZZ Euler(long r){
    long eu = to_long(1), p;

    for(p = 2; p * p <= r; p += 2){
        if(r % p == 0){
            eu *= p - 1;
            r /= p;

            while(r % p == 0){
                eu *= p;
                r /= p;
            }
        }

        if(p == 2)
        --p;
    }

    ZZ eu1 = to_ZZ(eu);

    // now r is prime or 1
    if(r == 1){
        return eu1;
    }
    else{
        return eu1 * (r - 1) ;
    }
}