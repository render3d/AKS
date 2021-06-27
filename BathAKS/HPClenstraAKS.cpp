/*
    A C++ implementation of Lenstra's algorithm, adapted from that from Hua Li, that uses the NTL library's faster arithemetic methods and basic thread pool functionality:
        https://researchportal.bath.ac.uk/en/publications/the-analysis-and-implementation-of-the-aks-algorithm-and-its-impr
*/

#include <math.h> // standard libraries
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
// #include <windows.h>
#include <unistd.h>
#include <stdio.h>
#include <cstdio>
#include <stdlib.h>
#include <cstdlib>
// #include <mmsystem.h>
#include <time.h>
#include <ctime>
#include <chrono>
#include <string>
#include <thread>
#include <sys/time.h>
#include <sys/resource.h>

#include "NTL/ZZ.h" // NTL Libraries
#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
#include "NTL/ZZX.h"
#include "NTL/vec_ZZ.h"
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/BasicThreadPool.h>
NTL_CLIENT

#include "PerfectPower.h" //Each Indepedent Test
#include "Euler.h"
#include "Congruence.h"

std::string getTime() {
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d-%H-%M-%S");
    auto time = oss.str();

    // std::cout << time;

    return time;
}

std::string getFilename() {
    std::string prfx = "log-LenstraHPC-";
    std::string sffx = getTime();
    std::string extn = ".txt";

    std::string filename = prfx + sffx + extn;

    return filename;
}

const auto ncores = std::thread::hardware_concurrency(); // machine cores - may return 0 when not able to detect
const auto SetNumThreads(ncores); // number of threads - should correspond to the number of available cores on your machine

std::string filename = getFilename();
std::ofstream perflog(filename, std::ios::app); // output result into file

bool lenstraAKS () {

    start:
    ZZ n;
    n = 0;

    std::cout << "Enter a positive integer number n you want to be tested:\n";
    std::cin >> n;

    if(n < 1){
        std::cout << "Integer n needs to be positive.\n\n";
        goto start;
    }
    else if(n == 1){
        std::cout << "1 is neither prime or composite.\n\n";
        return false;
    }
    else if(n == 2){
        std::cout << "2 is prime.\n\n";
        return true;
    }
    else if(n == 3){
        std::cout << "3 is prime.\n\n";
        return true;
    }

    std::cout << "n = " << n << "\n";
    perflog << "n = " << n << "\n";

    // start timing
    auto start = std::chrono::steady_clock::now();
    std::cout << "CPU clocks per second = " << CLOCKS_PER_SEC * 1000 << "\n";

    // Test if n is a perfect power
    int PP = PerfectPower(n);

    // returns 1 if n is a perfect power, 0 otherwise;
    if(PP == 1){
        auto finish = std::chrono::steady_clock::now();
        auto duration = finish - start;

        std::cout << n << " is a perfect power, hence is not prime.\n\n";
        std::cout << "Time Taken:" << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() << " milliseconds\n\n";

        perflog << n << " is a perfect power, hence is not prime.\n\n";
        perflog << "Time Taken:" << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() << " milliseconds\n\n";
        // perflog.close();
        return false;
    }

    // Find a suitable r
    ZZ r = to_ZZ(2);
    ZZ R;
    ZZ r1;

    while(r < n){ // line 3 of Fig 2.2
        ZZ R = GCD(r, n);
        if(R != 1 ){ // line 4 of Fig 2.2
            auto finish = std::chrono::steady_clock::now();
            auto duration = finish - start;

            std::cout << n << " has factors other than n and 1, hence is composite.\n\n";
            std::cout << R << " is a divisor.\n\n";
            std::cout << "Time Taken:" << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() << " milliseconds\n\n";

            perflog << n << " is composite.\n";
            perflog << R << " is a divisor.\n\n";
            perflog << "Time Taken:" << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() << " milliseconds\n\n";
            // perflog.close();
            return false;
        }
        else { // lines 5 and 6 of Fig 2.2
            ZZ v = to_ZZ(floor(power_long(to_long(log(n)), 2)));

            // order of n mod r is bigger than v;
            int p = 0;
            ZZ_p::init(r); // calculate mod r

            while(v <= r){
                ZZ x = to_ZZ(power_long(to_long(n), to_long(v))); // calculates x = n^v
                ZZ_p z = to_ZZ_p(x);
                if(z == to_ZZ_p(1)){
                    r1 = r; // store value of r;
                    r = n + 1;
                    break;
                }
                else{
                    v = v + 1;
                }
            }
        }
        r = r + 1; // line 7 of Fig 2.2
    }

    r = r1;
    std::cout << "r = " << r << "\n";
    perflog << "r = " << r << "\n";

    // calculate lines 11-13 of Fig 2.2
    ZZ r2 = Euler(to_long(r));
    std::cout << "Euler(" << r << ") = " << r2 << "\n";
    perflog << "Euler(" << r << ") = " << r2 << "\n";

    NTL_EXEC_RANGE(to_long(r2 - 1), first, last);

        for(long a = first; a <= last; ++a){ // line 9 of Fig 2.2
            int f = Congruence (a, n, r); // line 10 of Fig 2.2, returns 1 if condition holds, 0 otherwise

            if(f == 0){
                auto finish = std::chrono::steady_clock::now();
                auto duration = finish - start;

                std::cout << n << " is not prime.\n\n";
                std::cout << "the a which fails is " << a << "\n";
                std::cout << "Time Taken:" << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() << " milliseconds\n\n";

                perflog << n << " is not prime.\n"; // line 12 fails for particular a
                perflog << "the a which fails is " << a << "\n";
                perflog << "Time Taken:" << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() << " milliseconds\n\n";
                // perflog.close();

                return false;
            }
        }

    NTL_EXEC_RANGE_END;

    auto finish = std::chrono::steady_clock::now();
    auto duration = finish - start;

    std::cout << n << " is prime.\n\n";
    std::cout << "Time Taken:" << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() << " milliseconds\n\n";

    //n must be prime if went through this stage, output result to file
    perflog << "n is prime.\n";
    perflog << "Time Taken:" << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() << " milliseconds\n\n";
    // perflog.close();

    return true;
}

int main (int argc, char * argv[]){

    while (1) {
        int opt; // Give user opportunity to continue or exit program

        std::cout << "Press '1' to test a number, '0' to exit the program:\n";
        std::cin >> opt;

        if (opt == 1) {
            bool PRIME = lenstraAKS();
        }
        else if (opt == 0) {
            return(0); // Exit program
        }
    }
}
