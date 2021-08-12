// type test for chrono library steady_clock
// cp test/chronoTest.cpp test/HPCchronoTest.cpp
// $ g++ -g -O2 -std=c++11 -pthread -march=native test/chronoTest.cpp -o test/chronoTest.out -lntl -lgmp -lm

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
NTL_CLIENT

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

int main() {
    std::cout << getTime() << "\n";
    std::cout << getDate();

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

    std::cout << "Stop\n";
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() << '\n';
}
