/*
    A C++ implementation of Lenstra's algorithm, adapted from that from Hua Li:
        https://researchportal.bath.ac.uk/en/publications/the-analysis-and-implementation-of-the-aks-algorithm-and-its-impr
    Compile with:
        $ g++ -g -O2 -std=c++11 -pthread -march=native dir/foo.cpp -o dir/foo.out -lntl -lgmp -lm
        $ nvcc -o devVR/LenstraZnx.out devVR/LenstraZnx.cu
        $ nvcc devVR/LenstraZnx.cu -o devVR/LenstraZnx.out
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
#include <filesystem>
// #include <mmsystem.h>
#include <time.h>
#include <ctime>
#include <chrono>
#include <string>
#include <thread>
#include <array>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/types.h>

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

#include "PerfectPower.h" //Each Independent Test
#include "Euler.h"

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

std::string getDateTime() {
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d-%H-%M-%S");
    auto datetime = oss.str();

    return datetime;
}

std::string getFilename() {
    // std::string make = "/logs/" + getDate();
    // int result = mkdir(make.c_str(), 0777);
    std::string fldr = "logs/HPCznx/";
    // std::string fldr = "logs/" + getDate() + "/";
    std::string prfx = "LnstrZnx-";
    std::string sffx = getDateTime();
    std::string extn = ".csv";

    std::string filename = fldr + prfx + sffx + extn;

    return filename;
}

// std::filesystem::create_directory("logs/" + getDate());
unsigned int ncores = std::thread::hardware_concurrency(); // machine cores - may return 0 when not able to detect
const auto SetNumThreads(ncores); // number of threads - should correspond to the number of available cores on your machine

std::string filename = getFilename();
std::ofstream perflog(filename, std::ios::app); // output result into file

inline void fileWrite(const ZZ& n, const unsigned int& cores, const bool& PRIME, const long& time, const std::string& other) {
    perflog << n << "," << cores << "," << PRIME  << "," << time << "," << other << "\n";
}

__global__ void CongruenceZnx (ZZ *n, ZZ *r, ZZ *r2, long *u, long a) { // __global__ void kernel(ZZ *d_n, ZZ *d_r, ZZ *d_r2){
    // congruence test of polynomials in regular form

    // Thread indexing
    long i = to_long(threadIdx.x);

    // Perform this operation for every thread
    if (i < a) {                            // ensures kernel does not execute more threads than size of a
        ZZ_p::init(n);                      // initialise mod n

        ZZ_pX b = ZZ_pX(to_long(r), 1) - 1; // b = x^r - 1 (mod n);
        ZZ_pX e = ZZ_pX(1, 1);              // e = x (mod n)
        ZZ_pX d = PowerMod(e, n, b);        // d = x^n (mod b, n)

        ZZ_pX c = ZZ_pX(1, 1) - i;          // c = x - a (mod n);
        ZZ_pX f = PowerMod(c, n, b);        // f = (x - a)^n (mod b, n) - LHS
        ZZ_pX g = d - i;                    // g = x^n - a (mod b, n) - RHS

        if(f != g){
            u[i] = 0; // n is not prime
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

    std::cout << "n = " << n << "\n\n";

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

    h_r = r1;
    std::printf("r = %ld\n",to_long(h_r));

    ZZ h_r2 = Euler(to_long(h_r));
    std::printf("Euler(%ld) = %ld\n",to_long(h_r),to_long(h_r2));

    // Declare variables
    long h_an;
    long h_av[a];
    ZZ h_n = n;

    long* d_an;
    long* d_av;
    ZZ* d_n;
    ZZ* d_r;
    ZZ* d_r2;

    // Initialise variables
    h_an = to_long(h_r2 - 1);
    for (long x = 0; x < a; ++x) {
        h_av[x] = x + 1;
    }

    // Allocate memory on the device -- cudaMalloc(Location of Memory on Device,sizeof(int));
    cudaMalloc((void**)&d_an,sizeof(long));
    cudaMalloc((void**)&d_av,h_an*sizeof(long));
    cudaMalloc((void**)&d_n,sizeof(ZZ));
    cudaMalloc((void**)&d_r,sizeof(ZZ));
    cudaMalloc((void**)&d_r2,sizeof(ZZ));

    // Copy data from Host to Device
    cudaMemcpy(d_an,h_an,sizeof(long),cudaMemcpyHostToDevice);
    cudaMemcpy(d_av,h_av,h_an*sizeof(long),cudaMemcpyHostToDevice);
    cudaMemcpy(d_n,h_n,sizeof(ZZ),cudaMemcpyHostToDevice);
    cudaMemcpy(d_r,h_r,sizeof(ZZ),cudaMemcpyHostToDevice);
    cudaMemcpy(d_r2,h_r2,sizeof(ZZ),cudaMemcpyHostToDevice);

    // Configuration Parameters
    dim3 grid_size(1);
    dim3 block_size(h_an); // a threads in block

    // Launch Kernel -- CongruenceZnx<<<grid_size,block_size>>>(d_n,d_r,d_r2,d_av,d_an)
    CongruenceZnx<<<grid_size,block_size>>>(d_n,d_r,d_r2,d_av,d_an);

    // Copy data back to host
    cudaMemcpy(h_av,d_av,h_an*sizeof(long),cudaMemcpyDeviceToHost);

    // De-allocate memory
    cudaFree(d_an);
    cudaFree(d_av);
    cudaFree(d_n);
    cudaFree(d_r);
    cudaFree(d_r2);
    // free(h_av);

    for (long x = 0; x < h_an; ++x) {
        if(h_av[x] == 0){
            auto finish = std::chrono::steady_clock::now();
            auto duration = finish - start;
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();

            long a = x + 1;
            std::printf("%ld is not prime.\n",to_long(n));
            std::printf("The a which fails is %ld\n",a);
            std::printf("Time taken: %ld milliseconds\n",time);

            std::string note = "a = " + std::to_string(a) + "; r = " + std::to_string(to_long(r)) + "; phi(r) = " + std::to_string(to_long(r2));
            fileWrite(n,ncores,false,time,note);

            return false;
            // break;
        }
    }

    auto finish = std::chrono::steady_clock::now();
    auto duration = finish - start;
    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();

    std::printf("%ld is prime.\n",to_long(n));
    std::printf("Time taken: %ld milliseconds\n",time);

    std::string note = "n/a";
    fileWrite(n,ncores,true,time,note);

    return true;
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