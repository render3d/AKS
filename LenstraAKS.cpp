// Lenstra AKS implementation
// unoptimised program with no bignum dependencies
// Optimisation follows the completion of this handwritten implementation

//  if n is a prime power then return NOT PRIME

//  Find an integer r for which the multiplicative order of n (mod r) > round(lg^2 n)

/*  for each a ∈ [1, φ(r) − 1] do
        if (x + a)^n != x^n + a (mod n, x^r − 1) then return NOT PRIME
    return PRIME
*/