// calculates line 10 of Fig 2.2

int Congruence(long a, ZZ n, ZZ r){
    ZZ p::init(n); //mod n
    ZZ pX b = ZZ pX(to_long(r), 1) - 1; // b = x^r - 1;
    ZZ pX c = ZZ pX(1, 1) - a ; // c = x - a;
    ZZ pX f = PowerMod(c, n, b); // f =(x-a )^n mod c, n which is the RHS
    ZZ pX e = ZZ pX(1, 1);
    ZZ pX g = PowerMod(e, n, b); // x^n mod b, n
    g = g - a ; // g1 = x ^n - a mod c, n.

    if(f == g){
        return(1); // n is prime
    }
    else{
        return(0); // n is not prime.
    }
}