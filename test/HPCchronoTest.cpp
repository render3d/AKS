// type test for chrono library steady_clock
// cp test/chronoTest.cpp test/HPCchronoTest.cpp

#include <iostream>
#include <iomanip>
#include <ctime>
#include <chrono>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <typeinfo>
#include <string>

template <typename T> std::string type_name();

std::string getTime() {
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    std::ostringstream oss;
    oss << std::put_time(&tm, "%d-%m-%Y-%H-%M-%S");
    auto time = oss.str();

    std::cout << time << std::endl;

    return time;
}

int main() {
    getTime();

    auto now = std::chrono::high_resolution_clock::now();
    // auto now = std::chrono::steady_clock::now();
    std::cout << "Start\n";

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
