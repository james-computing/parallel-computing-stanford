#include <iostream>
#include <cmath>

#include "../include/2-lec.hpp"

int main() {
    const int N {10};
    const int terms {10};
    float * x = new float[N];
    float * y = new float[N];

    for (std::size_t i {0}; i < N; i++) {
        x[i] = 0.1* i;
    }

    sinx(N, terms, x, y);

    for (std::size_t i {0}; i < N; i++) {
        std::cout << y[i] - std::sin(x[i]) << std::endl;
    }

    delete[] x;
    delete[] y;
    x = nullptr;
    y = nullptr;

    return 0;
}