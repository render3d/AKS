
#include <cpuid.h>
#include <iostream>
#include <map>
#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <thread>
#include <filesystem>
#include <experimental/filesystem>

using namespace std;

//may return 0 when not able to detect
const auto nthreads = std::thread::hardware_concurrency();
char tmp[256];

struct CPUVendorID {
    unsigned int ebx;
    unsigned int edx;
    unsigned int ecx;

    string toString() const {
        return string(reinterpret_cast<const char *>(this), 12);
    }
};

int main() {
    unsigned int level = 0;
    unsigned int eax = 0;
    unsigned int ebx;
    unsigned int ecx;
    unsigned int edx;

    __get_cpuid(level, &eax, &ebx, &ecx, &edx);

    CPUVendorID vendorID { .ebx = ebx, .edx = edx, .ecx = ecx };

    map<string, string> vendorIdToName;
    vendorIdToName["GenuineIntel"] = "Intel";
    vendorIdToName["AuthenticAMD"] = "AMD";
    vendorIdToName["CyrixInstead"] = "Cyrix";
    vendorIdToName["CentaurHauls"] = "Centaur";
    vendorIdToName["SiS SiS SiS "] = "SiS";
    vendorIdToName["NexGenDriven"] = "NexGen";
    vendorIdToName["GenuineTMx86"] = "Transmeta";
    vendorIdToName["RiseRiseRise"] = "Rise";
    vendorIdToName["UMC UMC UMC "] = "UMC";
    vendorIdToName["Geode by NSC"] = "National Semiconductor";

    string vendorIDString = vendorID.toString();

    auto it = vendorIdToName.find(vendorIDString);
    string vendorName = (it == vendorIdToName.end()) ? "Unknown" : it->second;

    cout << "Max instruction ID: " << eax << endl;
    cout << "Vendor ID: " << vendorIDString << endl;
    cout << "Vendor name: " << vendorName << endl;
    cout << "Thread count: " << nthreads << endl;
    cout << "Current Path:"  << getcwd(tmp, 256) << endl;
}