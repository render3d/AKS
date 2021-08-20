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

#include "../devVR/CongruenceZ.h" // personal header

template <typename T> std::string type_name();

std::string getDate() {
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d");
    auto date = oss.str();

    return date;
}

std::string getTime() {
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    std::ostringstream oss;
    oss << std::put_time(&tm, "%H-%M-%S");
    auto time = oss.str();

    return time;
}

int main(int argc, char * argv[]) {

    ZZX f;

    std::cout << "Enter a polynomial using its coefficients in the place of the numbers";
    std::cout << " in the list [0,1,2,...,n] where n is its exponent (e.g. [2 -1 1] ";
    std::cout << " = x^2 - x + 2):\n";
    cin >> f;
    // test cases
    // [2 10 14 6] = 2 + 10x + 14x^2 + 6x^3
    // [9 16 7 1]  = 9 + 16x + 07x^2 + 1x^3
    // [14 10 6 2]
    // [16 9 7 1]
    // [2 -1 1]
    // [2 -1 -1]

    // printf("\nStart\n");

    // // *** DO STUFF HERE ***

    // ZZ coeffMaxF = getMaxCoeff(f);
    // std::cout << "\nMax coeff of f(x) = " << coeffMaxF << "\n";

    // ZZ x = ZZ(6);
    // ZZ resEval = evaluate(f,x);
    // std::cout << "\nWhen x = " << x << ", f(x) = " << resEval << "\n";

    ZZX g;
    std::cout << "Enter the polynomial you wish to multiply it by:\n";
    cin >> g;

    // auto then = std::chrono::high_resolution_clock::now();
    double then = GetTime();

    ZZX fgNTL = f * g;
    double lap1 = GetTime();
    // auto lap1 = std::chrono::high_resolution_clock::now();

    ZZX fgBin = polyMultiply(f,g);
    double lap2 = GetTime();
    // auto lap2 = std::chrono::high_resolution_clock::now();

    // auto tNTL = std::chrono::duration_cast<std::chrono::milliseconds>(lap1-then).count();
    double tNTL = lap1-then;
    std::cout << "Binary Segmentation: f(x) * g(x) = " << fgBin << " (" << tNTL*1000 << " milliseconds)\n";
    // auto tBin = std::chrono::duration_cast<std::chrono::milliseconds>(lap2-lap1).count();
    double tBin = lap2-lap1;
    std::cout << "NTL Polynomial Prod: f(x) * g(x) = " << fgNTL << " (" << tBin*1000 << " milliseconds)\n";

    if (fgNTL == fgBin) {
        printf("\nBinary Segmentation Multiply Successful.\n");
    }
    else {
        printf("\nBinary Segmentation Multiply Failed.\n");
    }
}
