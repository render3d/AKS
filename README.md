An implementation of the Lenstra variant of the AKS primailty test and a further, unimplemented variant suggestion from Crandall & Papadopoulos (2003): https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.95.58

Stages of Development:
1. *Basic implementation with no bignum dependencies - LenstraAKS.cpp
2. Optimised implementation, using MIRACL or NTL
3. GPU accelerated version in CUDA C++ for deployment on a GEFORCE RTX 3090 GPU cloud

Personal System Specifications for home benchmarking purposes:
- Eight-core AMD Ryzen 7 4700U, 8GB RAM
- Dual-core Intel(R) Pentium(R) CPU G860 @ 3.00GHz, 16GB DDR3 RAM @ 1333MHz, Radeon HD 7770 GPU