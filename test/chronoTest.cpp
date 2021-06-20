// type test for chrono library steady_clock

#include <iostream>
#include <chrono>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <typeinfo>

template <typename T> std::string type_name();

int main() {
    auto now = std::chrono::high_resolution_clock::now();
    std::cout << "Start\n";

    //do stuff here
    unsigned int microsecond = 1000000;
    usleep(3 * microsecond);//sleeps for 3 second
    std::cout << typeid(now).name() << '\n';

    auto then = std::chrono::high_resolution_clock::now();
    std::cout << "Stop\n";
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(then-now).count() << '\n';
}