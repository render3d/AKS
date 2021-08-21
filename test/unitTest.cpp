// type test for chrono library steady_clock
// cp test/chronoTest.cpp test/HPCchronoTest.cpp
// $ g++ -g -O2 -std=c++11 -pthread -march=native test/unitTest.cpp -o test/unitTest.out -lntl -lgmp -lm

#include <typeinfo>

#include <math.h> // standard libraries
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
// #include <windows.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
// #include <mmsystem.h>
#include <time.h>
#include <ctime>
#include <chrono>
#include <string>

#include "NTL/ZZ.h" // NTL Libraries
#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
#include "NTL/ZZX.h"
#include "NTL/vec_ZZ.h"
#include "NTL/ZZXFactoring.h"
NTL_CLIENT

#include "../devVR/CongruenceZ.h" // personal headers
#include "../devVR/CongruenceZp.h"

int main(int argc, char * argv[]) {

    // test cases
    // [2 10 14 6] = 2 + 10x + 14x^2 + 6x^3
    // [9 16 7 1]  = 9 + 16x + 07x^2 + 1x^3
    // [14 10 6 2]
    // [16 9 7 1]
    // [2 -1 1]
    // [2 -1 -1]

    printf("\nStart\n");

    // // *** DO STUFF HERE ***

    // ############### Multiply Tests ############################################

    // ZZX f;
    // std::cout << "Enter a polynomial using its coefficients in the place of the numbers";
    // std::cout << " in the list [0,1,2,...,n] where n is its exponent (e.g. [2 -1 1] ";
    // std::cout << " = x^2 - x + 2):\n";
    // cin >> f;

    // ZZX g;
    // std::cout << "Enter the polynomial you wish to multiply it by:\n";
    // cin >> g;

    // ZZX h;
    // std::cout << "And the polynomial you wish to mod the product by:\n";
    // cin >> h;

    // double then = GetTime();

    // ZZX fgNTL = f * g;
    // double lap1 = GetTime();

    // ZZX fgBin = polyMultiply(f,g);
    // double lap2 = GetTime();

    // double tNTL = lap1-then;
    // std::cout << "NTL Polynomial Prod: f(x) * g(x) = " << fgNTL << " (" << tNTL*1000 << " milliseconds)\n";

    // double tBin = lap2-lap1;
    // std::cout << "Binary Segmentation: f(x) * g(x) = " << fgBin << " (" << tBin*1000 << " milliseconds)\n";

    // if (fgNTL == fgBin) {
    //     printf("\nBinary Segmentation Multiply Successful.\n");
    // }
    // else {
    //     printf("\nBinary Segmentation Multiply Failed.\n");
    // }

    // ############### MulMod Tests ############################################

    // ZZX f;
    // std::cout << "Enter a polynomial using its coefficients in the place of the numbers";
    // std::cout << " in the list [0,1,2,...,n] where n is its exponent (e.g. [2 -1 1] ";
    // std::cout << " = x^2 - x + 2):\n";
    // cin >> f;

    // ZZX g;
    // std::cout << "Enter the polynomial you wish to multiply it by:\n";
    // cin >> g;

    // ZZX h;
    // std::cout << "And the polynomial you wish to mod the product by:\n";
    // cin >> h;

    // double then = GetTime();

    // ZZX fgnNTL = MulMod(f,g,h);
    // double lap1 = GetTime();

    // ZZX fgnBin = polyMulMod(f,g,h);
    // double lap2 = GetTime();

    // double tNTL = lap1-then;
    // std::cout << "NTL Polynomial Prod: f(x) * g(x) (mod h(x)) = " << fgnNTL << " (" << tNTL*1000 << " milliseconds)\n";

    // double tBin = lap2-lap1;
    // std::cout << "Binary Segmentation: f(x) * g(x) (mod h(x)) = " << fgnBin << " (" << tBin*1000 << " milliseconds)\n";

    // if (fgnNTL == fgnBin) {
    //     printf("\nBinary Segmentation MulMod Successful.\n");
    // }
    // else {
    //     printf("\nBinary Segmentation MulMod Failed.\n");
    // }

    // ############### PowMod Tests ############################################

    // ZZ_p::init(ZZ(2));

    // ZZX f;
    // std::cout << "Enter a polynomial using its coefficients in the place of the numbers";
    // std::cout << " in the list [0,1,2,...,n] where n is its exponent (e.g. [2 -1 1] ";
    // std::cout << " = x^2 - x + 2):\n";
    // cin >> f;

    // ZZ e;
    // std::cout << "Enter the exponent you wish to raise it to:\n";
    // cin >> e;

    // ZZX g;
    // std::cout << "And the polynomial you wish to mod the product by:\n";
    // cin >> g;

    // double then = GetTime();
    // ZZX fegBin = polyPowMod(f,e,g);
    // double lap1 = GetTime();

    // ZZ_pX f_p = to_ZZ_pX(f);
    // ZZ_pX g_p = to_ZZ_pX(g);

    // double lap2 = GetTime();
    // ZZ_pX fegNTL = PowerMod(f_p,e,g_p);
    // double lap3 = GetTime();

    // double tNTL = lap3-lap2;
    // std::cout << "NTL Polynomial PowerMod: f(x)^e (mod h(x)) = " << fegNTL << " (" << tNTL*1000 << " milliseconds)\n";

    // double tBin = lap1-then;
    // std::cout << "Binary Sgmtn polyPowMod: f(x)^e (mod h(x)) = " << fegBin << " (" << tBin*1000 << " milliseconds)\n";

    // if (fgnNTL == fgnBin) {
    //     printf("\nBinary Segmentation PowMod Successful.\n");
    // }
    // else {
    //     printf("\nBinary Segmentation PowMod Failed.\n");
    // }

    // ############### PowMod ZZ_p Tests ############################################

    ZZ_p::init(ZZ(2));

    ZZ_pX f;
    std::cout << "Enter a polynomial using its coefficients in the place of the numbers";
    std::cout << " in the list [0,1,2,...,n] where n is its exponent (e.g. [2 -1 1] ";
    std::cout << " = x^2 - x + 2):\n";
    cin >> f;

    ZZ e;
    std::cout << "Enter the exponent you wish to raise it to:\n";
    cin >> e;

    ZZ_pX g;
    std::cout << "And the polynomial you wish to mod the product by:\n";
    cin >> g;

    std::cout << "\nOperation " << f << "^" << e << " mod " << g << "\n";

    double then = GetTime();
    ZZ_pX fegBin = ZZpPowMod(f,e,g);
    double lap1 = GetTime();

    double lap2 = GetTime();
    ZZ_pX fegNTL = PowerMod(f,e,g);
    double lap3 = GetTime();

    double tNTL = lap3-lap2;
    std::cout << "NTL Polynomial PowerMod: f(x)^e (mod h(x)) = " << fegNTL << " (" << tNTL*1000 << " milliseconds)\n";

    double tBin = lap1-then;
    std::cout << "Binary Sgmtn polyPowMod: f(x)^e (mod h(x)) = " << fegBin << " (" << tBin*1000 << " milliseconds)\n";

    if (fegNTL == fegBin) {
        printf("\nBinary Segmentation PowMod Successful.\n");
    }
    else {
        printf("\nBinary Segmentation PowMod Failed.\n");
    }

    // test cases
    // [2 10 14 6] = 2 + 10x + 14x^2 + 6x^3
    // [9 16 7 1]  = 9 + 16x + 07x^2 + 1x^3
    // [14 10 6 2]
    // [16 9 7 1]
    // [2 -1 1]
    // [2 -1 -1]

    // [1 0 0 1]
    // [-1 1 1 0 -1 0 1]

}
