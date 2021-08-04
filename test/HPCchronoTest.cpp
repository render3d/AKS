// type test for chrono library steady_clock
// cp test/chronoTest.cpp test/HPCchronoTest.cpp
// g++ -g -O2 -std=c++11 -pthread -march=native test/HPCchronoTest.cpp -o test/test.out -lntl -lgmp -lm
// g++ test/HPCchronoTest.cpp -o test/HPCchronoTest.out

#include <typeinfo>
#include <string_view>

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

string getTime() {
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d-%H-%M-%S");
    auto time = oss.str();

    // std::cout << time;

    return time;
}

// std::ofstream perflog("logTest.txt", std::ios::app); // output result to txt file
std::ofstream perflog("logTest.csv", std::ios::app); // output result to csv file

inline void fileWrite(const int& n, const unsigned int& cores, const bool& PRIME, const long& time, const std::string& other) {
    perflog << n << "," << cores << "," << PRIME  << "," << time << "," << other << "\n";
}

int main () {
    std::cout << getTime();

    perflog << "Int, Cores, Prime (T/F), Time (milliseconds), Comments\n";

    auto now = std::chrono::high_resolution_clock::now();
    // auto now = std::chrono::steady_clock::now();
    std::cout << "\nStart\n";

    //do stuff here
    unsigned int microsecond = 1000000;
    usleep(3 * microsecond);//sleeps for 3 second
    // std::cout << typeid(now).name() << '\n';

    auto then = std::chrono::high_resolution_clock::now();
    // auto then = std::chrono::steady_clock::now();
    auto duration = then - now;
    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
    long l;

    std::cout << "Stop\n";
    std::cout << time << '\n';
    std::cout << typeid(time).name() << '\n';
    std::cout << typeid(l).name() << '\n';

    fileWrite(347,8,true,time,"note");
}
