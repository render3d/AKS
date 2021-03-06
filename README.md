# **University of Bath Computer Science MSc. Dissertation Project**
## **On the implementation of the AKS algorithm in CUDA**

An implementation of the Lenstra variant of the AKS primality test and a further, unimplemented variant suggestion from [Crandall & Papadopoulos (2003)](https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.95.58).

### Stages of Development:
- [ ] Review and validate (respectively) previous implementations from:
  - [ ]  [Gallot](http://yves.gallot.pagesperso-orange.fr/src/)
  - [ ]  [Hua Li](https://researchportal.bath.ac.uk/en/publications/the-analysis-and-implementation-of-the-aks-algorithm-and-its-impr) - Lenstra/$Z_n(x)$ version
- [ ] Implement Lenstra/$Z_n(x)$ version as GPU accelerated versions in CUDA C++ for deployment on a GEFORCE RTX 3090 GPU cloud.
- [ ] Develop $Z_n$ version from [Crandall & Papadopoulos (2003)](https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.95.58)
- [ ] Optimise above implementation, using MIRACL or NTL

### Other information
Personal System Specifications for home benchmarking purposes:
- Eight-core AMD Ryzen 7 4700U, 8GB RAM
- Dual-core Intel(R) Pentium(R) CPU G860 @ 3.00GHz, 16GB DDR3 RAM @ 1333MHz, Radeon HD 7770 GPU