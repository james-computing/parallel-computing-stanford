#include "../include/2-lec.hpp"

#include <iostream>

// Computes sin(value) for each value in an array x.
// It uses the Taylor expansion sin(x) = x -x^3/3! +x^5/5! ...
// This function is written in a way that is harder for hardware to parallelize.
// N = length of x
// terms = number of terms from the Taylor expansion
// x = array of inputs
// y = array of outputs
void sinx(size_t const N, size_t const terms, float const * const x, float * const y) {
    for (size_t i {0}; i < N; ++i) {
        float value = x[i];

        float numer = x[i] * x[i] * x[i];
        int denom  = 6; // 3!
        int sign = -1;

        for (size_t j {1}; j <= terms; ++j) {
            value += sign * numer / denom;
            numer *= x[i] * x[i];
            denom *= (2*j + 2) * (2*j + 3);
            sign *= -1;
        }

        y[i] = value;
    }
}

void my_thread_func(my_args const * const args) {
    sinx(args->N, args->terms, args->x, args->y);
}

// This computes the first half in the main thread and the second half in a new thread
void parallel_sinx(size_t const N, size_t const terms, float const * const x, float * const y) {
    my_args const args = my_args {
        .N {N/2},
        .terms {terms},
        .x {x},
        .y {y}
    };

    std::thread my_thread = std::thread(my_thread_func, &args); // lauch thread
    sinx(N-args.N, terms, x + args.N, y + args.N); // do work on main thread
    my_thread.join();
}