#pragma once

#include <cstddef> // size_t
#include <thread>
#include <immintrin.h> // SIMD
#include <cmath> // sine
#include <chrono> // measure time
#include <cstring> // for memcpy
#include <iostream>

namespace lec2
{

void sinx(size_t const N, size_t const terms, float const * const x, float * const y);

struct my_args {
    size_t const N;
    size_t const terms;
    float const * const x;
    float * const y;
};

void my_thread_func(my_args const * const args);
void parallel_sinx(size_t const N, size_t const terms, float const * const x, float * const y);

void simd_example();
#ifdef AVX512
void simd_with_branching_AVX512(size_t const N, float const * const x, float * const y);
#endif
void bad_init_for_simd(size_t const N, float * & x);
void good_init_for_simd(size_t const N, float * & x);
void no_simd_baseline(size_t const N, float const * const x, float * const y);
void simd_with_branching_AVX(size_t const N, float const * const x, float * const y);
void simd_optimized_for_good_init(size_t const N, float const * const x, float * const y);
bool are_close(size_t const N, float const * const y1, float const * const y2, float const tolerance);
void execute_examples();

// templates must be put in header file
template <typename T>
void print_array(size_t const N, T const * const a) {
    T const * p = a;
    for (size_t i {0}; i < N; ++i) {
        std::cout << *p << ' ';
        ++p;
    }
    std::cout << std::endl;
}
}