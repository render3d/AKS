/*
    A C++ implementation of Lenstra's algorithm, adapted from that from Hua Li:
        https://researchportal.bath.ac.uk/en/publications/the-analysis-and-implementation-of-the-aks-algorithm-and-its-impr
    Compile with:
        $ g++ -g -O2 -std=c++11 -pthread -march=native dir/foo.cpp -o dir/foo.out -lntl -lgmp -lm
        $ nvcc -o HPCLenstraZnx.out HPCLenstraZnx.cu
        $ nvcc HPCLenstraZnx.cu -o HPCLenstraZnx.out
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

// includes CUDA
#include <cuda_runtime.h>

// includes, project
#include <helper_cuda.h>
#include <helper_functions.h> // helper functions for SDK examples

#include "PerfectPower.h" //Each Indepedent Test
#include "Euler.h"

std::string getTime() {
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d-%H-%M-%S");
    auto time = oss.str();

    return time;
}

std::string getFilename() {
    std::string prfx = "log-LenstraZ-";
    std::string sffx = getTime();
    std::string extn = ".csv";

    std::string filename = prfx + sffx + extn;

    return filename;
}

unsigned int ncores = std::thread::hardware_concurrency(); // machine cores - may return 0 when not able to detect
const auto SetNumThreads(ncores); // number of threads - should correspond to the number of available cores on your machine

std::string filename = getFilename();
std::ofstream perflog(filename, std::ios::app); // output result into file

inline void fileWrite(const ZZ& n, const unsigned int& cores, const bool& PRIME, const long& time, const std::string& other) {
    perflog << n << "," << cores << "," << PRIME  << "," << time << "," << other << "\n";
}

__global__ void CongruenceZnx (ZZ *n, ZZ *r, ZZ *r2) { // __global__ void kernel(ZZ *d_n, ZZ *d_r, ZZ *d_r2){
    // congruence test of polynomials in regular form

    // Thread indexing
    int i = threadIdx.x;

    // Perform this operation for every thread

    for(long a = 1; a <= to_long(*r2 - 1); ++a){
        ZZ_p::init(*n); //mod n
        ZZ_pX b = ZZ_pX(to_long(*r), 1) - 1; // b = x^r - 1;
        ZZ_pX c = ZZ_pX(1, 1) - a ; // c = x - a;
        ZZ_pX f = PowerMod(c, *n, b); // f =(x - a)^n mod c, n which is the RHS
        ZZ_pX e = ZZ_pX(1, 1);
        ZZ_pX g = PowerMod(e, *n, b); // x^n mod b, n
        g = g - a ; // g1 = x^n - a mod c, n.

        if(f == g){
            return(1); // n is prime
        }
        else{
            return(a); // n is not prime.
        }
    }
}

inline bool Lenstra (const ZZ& n) {
    if(n < 1){
        std::printf("Integer n needs to be positive.\n");
        return false;
    }
    else if(n == 1){
        std::printf("1 is neither prime or composite.\n");
        return false;
    }
    else if(n == 2){
        std::printf("2 is prime.\n");
        return true;
    }
    else if(n == 3){
        std::printf("3 is prime.\n");
        return true;
    }

    std::printf("n = %ld\n",to_long(n));

    // start timing
    auto start = std::chrono::steady_clock::now();

    // Test if n is a perfect power
    int PP = PerfectPower(n);

    // returns 1 if n is a perfect power, 0 otherwise;
    if(PP == 1){
        auto finish = std::chrono::steady_clock::now();
        auto duration = finish - start;
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();

        std::printf("%ld is not prime.\n",to_long(n));
        std::printf("%ld is a perfect power.\n",to_long(n));
        std::printf("Time taken: %ld milliseconds\n",time);

        std::string note = std::to_string(to_long(n)) + " is a perfect power";
        fileWrite(n,ncores,false,time,note);

        return false;
    }

    // Find a suitable r
    ZZ r = to_ZZ(2);
    ZZ R;
    ZZ r1;

    while(r < n){
        ZZ R = GCD(r, n);
        if(R != 1 ){
            auto finish = std::chrono::steady_clock::now();
            auto duration = finish - start;
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();

            std::printf("%ld is not prime.\n",to_long(n));
            std::printf("%ld is a divisor.\n",to_long(R));
            std::printf("Time taken: %ld milliseconds\n",time);

            std::string note = std::to_string(to_long(R)) + " is a divisor";
            fileWrite(n,ncores,false,time,note);

            return false;
            break;
        }
        else {
            ZZ v = to_ZZ(floor(power_long(to_long(log(n)), 2)));

            // order of n mod r is bigger than v;
            int p = 0;
            ZZ_p::init(r); // calculate mod r

            while(v <= r){
                ZZ x = to_ZZ(power_long(to_long(n), to_long(v))); // calculates x = n^v
                ZZ_p z = to_ZZ_p(x);
                if(z == to_ZZ_p(1)){
                    r1 = r; // store value of r
                    r = n + 1;
                    break;
                }
                else{
                    v = v + 1;
                }
            }
        }
        r = r + 1;
    }

    r = r1;
    std::printf("r = %ld\n",to_long(r));

    ZZ r2 = Euler(to_long(r));
    std::printf("Euler(%ld) = %ld\n",to_long(r),to_long(r2));

    // Declare variables
    int *h_f, *d_f;

    // Allocate memory on the device -- cudaMalloc(Location of Memory on Device,sizeof(int));
    cudaMalloc((void**)&d_f,sizeof(int));

    // Copy data from Host to Device
    cudaMemcpy(d_f,h_f,sizeof(int),cudaMemcpyHostToDevice);

    // Configuration Parameters
    dim3 grid_size(1);
    dim3 block_size(N); \\ N threads in block

    // Launch Kernel -- CongruenceZnx<<<grid_size,block_size>>>(d_n,d_r,d_r2)
    CongruenceZnx<<<grid_size,block_size>>>(n,r,r2);

    // Copy data back to host
    cudaMemcpy(h_f,d_f,sizeof(int),cudaMemcpyDeviceToHost);

    // De-allocate memory
    cudaFree(d_f);
    free(h_f);

    if(f == 1){
        auto finish = std::chrono::steady_clock::now();
        auto duration = finish - start;
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();

        std::printf("%ld is prime.\n",to_long(n));
        std::printf("Time taken: %ld milliseconds\n",time);

        std::string note = "n/a";
        fileWrite(n,ncores,true,time,note);

        return true;
    }
    else {
        auto finish = std::chrono::steady_clock::now();
        auto duration = finish - start;
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();

        std::printf("%ld is not prime.\n",to_long(n));
        std::printf("The a which fails is %ld\n",f);
        std::printf("Time taken: %ld milliseconds\n",time);

        std::string note = "a = " + std::to_string(f) + "; r = " + std::to_string(to_long(r)) + "; phi(r) = " + std::to_string(to_long(r2));
        fileWrite(n,ncores,false,time,note);

        return false;
        // break;
    }

}

int main (int argc, char * argv[]) {

    perflog << "Int, Cores, Prime (T/F), Time (milliseconds), Comments\n";

    bool prime;
    ZZ n;
    n = 0;

    std::printf("Enter a positive integer number n you want to be tested:\n");
    std::cin >> n;

    prime = Lenstra(n);

    return 0;
}