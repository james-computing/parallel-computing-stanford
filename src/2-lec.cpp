#include "../include/2-lec.hpp"

// Computes sin(value) for each value in an array x.
// It uses the Taylor expansion sin(x) = x -x^3/3! +x^5/5! ...
// This function is written in a way that is harder for hardware to parallelize.
// N = length of x
// terms = number of terms from the Taylor expansion
// x = array of inputs
// y = array of outputs
void sinx(const size_t N, const size_t terms, float * const x, float * const y) {
    for (size_t i {0}; i < N; i++) {
        float value = x[i];

        float numer = x[i] * x[i] * x[i];
        int denom  = 6; // 3!
        int sign = -1;

        for (size_t j {1}; j <= terms; j++) {
            value += sign * numer / denom;
            numer *= x[i] * x[i];
            denom *= (2*j + 2) * (2*j + 3);
            sign *= -1;
        }

        y[i] = value;
    }
}