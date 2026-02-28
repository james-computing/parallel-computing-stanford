#pragma once

#include <cstddef>
#include <thread>

void sinx(size_t const N, size_t const terms, float const * const x, float * const y);

struct my_args {
    size_t const N;
    size_t const terms;
    float const * const x;
    float * const y;
};

void my_thread_func(my_args const * const args);
void parallel_sinx(size_t const N, size_t const terms, float const * const x, float * const y);