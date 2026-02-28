#include <iostream>
//#include <cmath>
#include <chrono>

#include "../include/2-lec.hpp"

int main() {
    // Alloc resources
    size_t const N {100};
    size_t const terms {10};
    float * x = new float[N];
    float * y = new float[N];
    float * sin = new float[N];

    // Initialize x
    for (std::size_t i {0}; i < N; ++i) {
        x[i] = 0.1* (i % 10);
    }

    // Variables to measure time
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    // Compute sin(x[i]) for each i, and output the results to y
    begin = std::chrono::steady_clock::now();
    sinx(N, terms, x, y);
    end = std::chrono::steady_clock::now();
    
    std::cout << "time sinx: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count() << std::endl;

    // Compute again, but with 2 threads
    begin = std::chrono::steady_clock::now();
    parallel_sinx(N, terms, x, y);
    end = std::chrono::steady_clock::now();

    std::cout << "time parallel_sinx: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count() << std::endl;

    // Cleanup
    delete[] sin;
    delete[] y;
    delete[] x;
    sin = nullptr;
    y = nullptr;
    x = nullptr;

    return 0;
}