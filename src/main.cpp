#include <iostream>
#include <cmath>

#include "../include/2-lec.hpp"

int main() {
    // Alloc resources
    const int N {10};
    const int terms {10};
    float * x = new float[N];
    float * y = new float[N];

    // Initialize x
    for (std::size_t i {0}; i < N; i++) {
        x[i] = 0.1* i;
    }

    // Compute sin(x[i]) for each i, and output the results to y
    sinx(N, terms, x, y);

    // Compare y[i] with sin(x[i])
    for (std::size_t i {0}; i < N; i++) {
        std::cout << y[i] - std::sin(x[i]) << std::endl;
    }

    // Cleanup
    delete[] x;
    delete[] y;
    x = nullptr;
    y = nullptr;

    return 0;
}