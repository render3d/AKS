/*
    A C++ implementation of the original C++ algorithm by Hua Li that uses the NTL library;
        https;//researchportal.bath.ac.uk/en/publications/the-analysis-and-implementation-of-the-aks-algorithm-and-its-impr
*/

#include <math.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
// #include <windows.h>
// #include <mmsystem.h>
#include <time.h>
#include <chrono>

#include <NTL/ZZ.h> // NTL Libraries
#include <NTL/RR.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZX.h>
#include <NTL/vec_ZZ.h>

#include "PerfectPower.h" //Each Indepedent Test
#include "IsPrime.h"
#include "LargestPrime.h"
#include "Congruence.h"

int main(int argc, char argv[]){

    start:
    ZZ n;
    n = 0;

    std::cout << "Enter a positive integer number n you want to be tested \n";
    std::cin >> n;

    if(n < 1){
        std::cout << "Integer n needs to be positive \n\n";
        goto start;
    }
    else if(n == 1){
        std::cout << "1 is neither prime or composite . \n\n";
        goto start;
    }
    else if(n == 2){
        std::cout << "2 is prime. \n\n";
        goto start;
    }
    else if(n == 3){
        std::cout << "3 is prime. \n\n";
        goto start;
    }

    std::ofstream my_file ("AKS_original.doc", std::ios::app); // output result into file
    std::cout << "n = " << n << "\n";
    my_file << "n = " << n << "\n";

    // start timing
    auto start = std::chrono::steady_clock::now();
    std::cout << "CPU clocks per second = " << CLOCKS_PER_SEC * 1000 << "\n";

    int PP = PerfectPower(n); // Check if n is a perfect power; returns 1 if n is a perfect power, 0 otherwise

    if(PP==1){
        my_file << n << " is a perfect power, hence is not prime \n\n";
        auto finish = std::chrono::steady_clock::now();
        auto duration = finish - start;
        std::cout << "Time Taken = " << duration << " milliseconds \n\n";
        my_file << "Time Taken = " << duration << " milliseconds \n\n";
        my_file.close();
        int choice;

        // Give user opportunity to continue or exit program
        std::cout << "Press '1' to test a new number, '0' to exit the program\n";
        std::cin >> choice;

        if(choice == 1){
            goto start;
        }
        else if(choice == 0){
            return(0); // Exit program
        }
    }
    else{
        // continue
    }

    // Find a suitable r
    ZZ r = to_ZZ(2);

    while(r < n){ // line 3 of Fig 2.1
        ZZ r1 = GCD(r, n);

        if(r1 != 1){ // line 4 of Fig 2.1
            std::cout << "n has factors other than n and 1, hence is composite \n\n";

            my_file << n << " is composite. \n";
            my_file << r1 << " is a divisor \n\n";

            std::cout << r1 << " is a divisor \n\n";

            auto finish = std::chrono::steady_clock::now();
            auto duration = finish - start;

            std::cout << "Time Taken = " << duration << " milliseconds \n\n";

            my_file << "Time Taken = " << duration << " milliseconds \n\n";
            my_file.close();

            int choice;

            // Give user opportunity to continue or exit program
            std::cout << " Press '1' to test a new number, '0' to Exit the program\n"; std::cin >> choice;
            if(choice == 1){
                goto start;
            }
            else if(choice == 0){
                return(0); // Exit program
            }
        }

        int prime = IsPrime(r); // line 5 of Fig 2.1
        while(prime == 0){
            r = r + 1;// increase r until it becomes prime
            prime = IsPrime(r);
        }

        ZZ q = LargestPrime(r - 1); // line 6 of Fig 2.1, set q to be largest prime factor of r-1
        ZZ_p::init(r); // mod r
        double t1 = 4 * sqrt(to_long(r)) * (log(n)/log(2));

        // test on lest hand side in line 7 of Fig 2.1
        ZZ t2;
        t2 = power(n, to_long((r - 1) / q)) % r;

        // test on right hand side in line 7 of Fig 2.1
        if(to_double(q) >= t1 && t2 != 1){
            break; // line 8 of Fig 2.1, q satisfies required conditions, stop while loop
        }
        else{
            r = r + 1; // line 9 of Fig 2.1
        }
    }

    if(r > n){
        r = n;
    }

    std::cout << "r = " << r << "\n";
    my_file << "r = " << r << "\n";
    long a;
    long b = to_long(2 * sqrt(to_long(r)) * (log(n) / log(2)));

    for(a = 1; a <= b; a++){ // line 11 of Fig 2.1
        int f = Congruence(a, n, to_ZZ(r)); // line 12 of Fig 2.1

        // Returns 0 if condition fails, returns 1 if holds
        if(f == 0){
            my_file << "a fails at " << a << "\n";
            std::cout << "a fails at " << a << "\n";

            // Line 12 fails for an a, hence n is composite
            my_file << "Hence n is not prime \n";
            auto finish = std::chrono::steady_clock::now();
            auto duration = finish - start;

            std::cout << "Time Taken=" << duration << " milliseconds \n\n";
            my_file << "Time Taken=" << duration << " milliseconds \n\n";
            my_file.close();
            std::cout << n << " is not prime. \n\n";

            int choice; // Give user opportunity to continue or exit program
            std::cout << "Press '1' to test a new number, '0' to exit the program \n"; std::cin >> choice;
            if(choice == 1){
                goto start;
            }
            else if(choice == 0){
                return (0); // Exit program
            }
        }
    }

    auto finish = std::chrono::steady_clock::now();
    auto duration = finish - start;

    my_file << "n is prime \n"; // line 13 Fig 2.1
    std::cout << "Time Taken = " << duration << " milliseconds \n\n";
    my_file << "Time Taken = " << duration << " milliseconds \n\n";
    my_file.close();
    std::cout << n << " is prime. \n\n";

    int choice; // Give user opportunity to continue or Exit program
    std::cout << "Press '1' to test a new number, '0' to exit the program \n";
    std::cin >> choice;
    if(choice == 1){
        goto start;
    }
    else if(choice == 0){
        return(1); // Exit program
    }
}