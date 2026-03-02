#include "../include/lec2.hpp"

namespace lec2
{

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

//******************************************** SIMD ************************************

void simd_example() {
    // SIMD (single instruction, multiple data)
    // Make vectors of 8 floats to be used by SIMD
    __m256 const evens {_mm256_set_ps(2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0)};
    __m256 const odds {_mm256_set_ps(1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0)};
    // Subtract these vectors
    __m256 const result {_mm256_sub_ps(evens, odds)};
    // Cast to a pointer to access individual elements
    float const * const f  {(float const *) &result};

    // Print the results
    std::cout << "\nResult from simd_example:\n";
    for (char i{0}; i < 8; ++i) {
        std::cout << f[i] << " ";
    }
    std::cout << std::endl;
}

void bad_init_for_simd(size_t const N, float * & x) {
    // The worst possible case for the function simd_with_branching_AVX
    // is for the true condition too happen many times,
    // but being distributed with at most one case for each vector
    // of 8 floats.
    // This corresponds to a mask like [true, false, false, false, false, false, false, false].
    for (size_t i {0}; i < N; ++i) {
        if (i % 8 == 0) {
            x[i] = 1.0;
        } else {
            x[i] = -1.0;
        }
    }
}

void good_init_for_simd(size_t const N, float * & x) {
    // The bad case can be rearranged into a good case by grouping together
    // the positive values. The positive values take longer to execute,
    // by how simd_with_branching_AVX was designed.
    size_t i {0};
    for (; i < N/8; ++i) {
        x[i] = 1.0;
    }
    for (; i < N; ++i) {
        x[i] = -1.0;
    }
}

void no_simd_baseline(size_t const N, float const * const x, float * const y) {
    for (size_t i {0}; i < N; i++) {
        float t {x[i]};

        if (t > 0.0) {
            // Same from the lecture
            t = t*t;
            t = 50.0*t;
            t = t +100.0;
            // Some more operations to make it expensive
            t = t*t;
            t = t-30.0;
            t = 50.0*t;
            t = t+10.0;

            t = t*t;
            t = t-30.0;
            t = 50.0*t;
            t = t+10.0;

            t = t*t;
            t = t-30.0;
            t = 50.0*t;
            t = t+10.0;
        } else {
            // Same from the lecture
            t = t + 30.0;
            t = t/10.0;
        }

        y[i] = t;
    }
}

// Use if the cpu supports AVX512
#ifdef AVX512
// It assumes that N is multiple of 8
void simd_with_branching_AVX512(size_t const N, float const * const x, float * const y) {
    // Pointer to access the values of x and y, one at a time, using pointer arithmetic
    float const * input_ptr {x};
    float * output_ptr {y};

    // Vectors to be used by SIMD
    __m256 input_vector;
    __m256 result_vector;
    __mmask8 mask;

    // Pointer just for reading the results from SIMD
    float const * const result_float_ptr {(float const *) &result_vector};

    // Some constants to use in the SIMD operations
    __m256 const zero {_mm256_set1_ps(0.0)};
    __m256 const ten {_mm256_set1_ps(10.0)};
    __m256 const thirty {_mm256_set1_ps(30.0)};
    __m256 const fifty {_mm256_set1_ps(50.0)};
    __m256 const hundred {_mm256_set1_ps(100.0)};

    for (size_t i {0}; i < N/8; ++i) {
        // Copy 8 floats from x to v, starting from where input_ptr is pointing to.
        // Use the unaligned version, loadu.
        input_vector = _mm256_loadu_ps(input_ptr);

        // Compute the mask.
        // It is true when value > 0.0.
        
        mask = _mm256_cmp_ps_mask(input_vector, zero, _CMP_GT_OQ);
        
        // Do the operations corresponding to the condition value > 0.0
        // Do something expensive.
        // Same from lecture.   
        result_vector = _mm256_mask_mul_ps(input_vector, mask, input_vector, input_vector);
        result_vector = _mm256_mask_mul_ps(input_vector, mask, result_vector, fifty);
        result_vector = _mm256_mask_add_ps(input_vector, mask, result_vector, hundred);
        // More operations to make expensive
        result_vector = _mm256_mask_mul_ps(input_vector, mask, result_vector, result_vector);
        result_vector = _mm256_mask_sub_ps(input_vector, mask, result_vector, thirty);
        result_vector = _mm256_mask_mul_ps(input_vector, mask, result_vector, fifty);
        result_vector = _mm256_mask_add_ps(input_vector, mask, result_vector, ten);

        result_vector = _mm256_mask_mul_ps(input_vector, mask, result_vector, result_vector);
        result_vector = _mm256_mask_sub_ps(input_vector, mask, result_vector, thirty);
        result_vector = _mm256_mask_mul_ps(input_vector, mask, result_vector, fifty);
        result_vector = _mm256_mask_add_ps(input_vector, mask, result_vector, ten);
        
        result_vector = _mm256_mask_mul_ps(input_vector, mask, result_vector, result_vector);
        result_vector = _mm256_mask_sub_ps(input_vector, mask, result_vector, thirty);
        result_vector = _mm256_mask_mul_ps(input_vector, mask, result_vector, fifty);
        result_vector = _mm256_mask_add_ps(input_vector, mask, result_vector, ten);

        // Invert the mask
        mask = _knot_mask8(mask);

        // Do the operations corresponding to the condition value <= 0.0
        // Do something cheap.
        // Same from lecture.
        result_vector = _mm256_mask_add_ps(input_vector, mask, input_vector, thirty);
        result_vector = _mm256_mask_div_ps(input_vector, mask, input_vector, ten);

        // Copy the results to output
        memcpy(output_ptr, result_float_ptr, 8*sizeof(float));

        // Advance to next 8 float values
        input_ptr += 8;
        output_ptr += 8;
    }
}
#endif

// Do a version for the AVX family.
// In this case we don't have the mask.
// It assumes that N is multiple of 8
void simd_with_branching_AVX(size_t const N, float const * const x, float * const y) {
    // Pointer to access the values of x and y, one at a time, using pointer arithmetic
    float const * input_ptr {x};
    float * output_ptr {y};

    // Vectors to be used by SIMD.
    // In the absence of masks, use two result vectors, one for each condition.
    __m256 input_vector;
    __m256 result_true_vector;
    __m256 result_false_vector;

    // Pointers just for reading the results from SIMD
    float const * result_true_float_ptr {(float const *) &result_true_vector};
    float const * result_false_float_ptr {(float const *) &result_false_vector};

    // Some constants to use in the SIMD operations
    __m256 const ten {_mm256_set1_ps(10.0)};
    __m256 const thirty {_mm256_set1_ps(30.0)};
    __m256 const fifty {_mm256_set1_ps(50.0)};
    __m256 const hundred {_mm256_set1_ps(100.0)};

    for (size_t i {0}; i < N/8; ++i) {

        // Copy 8 floats from x to input_vector, starting from where input_ptr is pointing to.
        // Use the unaligned version, loadu.
        input_vector = _mm256_loadu_ps(input_ptr);

        // Do the operations corresponding to the condition value > 0.0
        // Do something expensive.
        // First use the same operations from the lecture.
        result_true_vector = _mm256_mul_ps(input_vector, input_vector);
        result_true_vector = _mm256_mul_ps(result_true_vector, fifty);
        result_true_vector = _mm256_add_ps(result_true_vector, hundred);
        // Then add some more operations to make it more expensive.
        result_true_vector = _mm256_mul_ps(result_true_vector, result_true_vector);
        result_true_vector = _mm256_sub_ps(result_true_vector, thirty);
        result_true_vector = _mm256_mul_ps(result_true_vector, fifty);
        result_true_vector = _mm256_add_ps(result_true_vector, ten);

        result_true_vector = _mm256_mul_ps(result_true_vector, result_true_vector);
        result_true_vector = _mm256_sub_ps(result_true_vector, thirty);
        result_true_vector = _mm256_mul_ps(result_true_vector, fifty);
        result_true_vector = _mm256_add_ps(result_true_vector, ten);
        
        result_true_vector = _mm256_mul_ps(result_true_vector, result_true_vector);
        result_true_vector = _mm256_sub_ps(result_true_vector, thirty);
        result_true_vector = _mm256_mul_ps(result_true_vector, fifty);
        result_true_vector = _mm256_add_ps(result_true_vector, ten);

        // Do the operations corresponding to the condition value <= 0.0
        // Do something cheap.
        // Use the same operations from the lecture.
        result_false_vector = _mm256_add_ps(input_vector, thirty);
        result_false_vector = _mm256_div_ps(result_false_vector, ten);

        // Copy the results to output
        for (char j {0}; j < 8; ++j) {
            if (*input_ptr > 0.0) {
                *output_ptr = *result_true_float_ptr;
            } else {
                *output_ptr = *result_false_float_ptr;
            }
            ++input_ptr;
            ++output_ptr;
            ++result_true_float_ptr;
            ++result_false_float_ptr;
        }

        // Reset pointers
        result_true_float_ptr = (float const *) &result_true_vector;
        result_false_float_ptr = (float const *) &result_false_vector;
    }
}

void simd_optimized_for_good_init(size_t const N, float const * const x, float * const y) {
    // Pointer to access the values of x and y, one at a time, using pointer arithmetic
    float const * input_ptr {x};
    float * output_ptr {y};

    // Vectors to be used by SIMD.
    __m256 input_vector;
    __m256 result_vector;

    // Pointers just for reading the results from SIMD
    float const * const result_float_ptr {(float const *) &result_vector};

    // Some constants to use in the SIMD operations
    __m256 const ten {_mm256_set1_ps(10.0)};
    __m256 const thirty {_mm256_set1_ps(30.0)};
    __m256 const fifty {_mm256_set1_ps(50.0)};
    __m256 const hundred {_mm256_set1_ps(100.0)};
    
    size_t i {0};
    for (; i < N/64; ++i) {
        // Copy 8 floats from x to input_vector, starting from where input_ptr is pointing to.
        // Use the unaligned version, loadu.
        input_vector = _mm256_loadu_ps(input_ptr);

        // Do the operations corresponding to the condition value > 0.0
        // Do something expensive.
        // First use the same operations from the lecture.
        result_vector = _mm256_mul_ps(input_vector, input_vector);
        result_vector = _mm256_mul_ps(result_vector, fifty);
        result_vector = _mm256_add_ps(result_vector, hundred);
        // Then add some more operations to make it more expensive.
        result_vector = _mm256_mul_ps(result_vector, result_vector);
        result_vector = _mm256_sub_ps(result_vector, thirty);
        result_vector = _mm256_mul_ps(result_vector, fifty);
        result_vector = _mm256_add_ps(result_vector, ten);

        result_vector = _mm256_mul_ps(result_vector, result_vector);
        result_vector = _mm256_sub_ps(result_vector, thirty);
        result_vector = _mm256_mul_ps(result_vector, fifty);
        result_vector = _mm256_add_ps(result_vector, ten);
        
        result_vector = _mm256_mul_ps(result_vector, result_vector);
        result_vector = _mm256_sub_ps(result_vector, thirty);
        result_vector = _mm256_mul_ps(result_vector, fifty);
        result_vector = _mm256_add_ps(result_vector, ten);

        // Copy the results to output
        memcpy(output_ptr, result_float_ptr, 8*sizeof(float));

        // Advance to next 8 values
        input_ptr += 8;
        output_ptr += 8;
    }

    for (; i < N/8; ++i) {
        // Copy 8 floats from x to input_vector, starting from where input_ptr is pointing to.
        // Use the unaligned version, loadu.
        input_vector = _mm256_loadu_ps(input_ptr);
        
        // Do the operations corresponding to the condition value <= 0.0
        // Do something cheap.
        // Use the same operations from the lecture.
        result_vector = _mm256_add_ps(input_vector, thirty);
        result_vector = _mm256_div_ps(result_vector, ten);

        // Copy the results to output
        memcpy(output_ptr, result_float_ptr, 8*sizeof(float));

        // Advance to next 8 values
        input_ptr += 8;
        output_ptr += 8;
    }
}

// Equality is obtained by taking tolerance = 0
bool are_close(size_t const N, float const * const y1, float const * const y2, float const tolerance) {
    for (size_t i {0}; i < N; ++i) {
        if (std::abs(y1[i] - y2[i]) > tolerance){
            std::cout << "\nFailed to be close: y1["<< i << "] = " << y1[i] << ", y2[" << i << "] = " << y2[i] << std::endl; 
            return false;
        }
    }
    return true;
}

void execute_examples() {
    // Alloc resources
    size_t const exponent {25};
    size_t const N {1<<exponent}; // 2^exponent
    size_t const terms {10};
    float const tolerance {1E-7}; // to compare with sin from cmath library
    float * x = new float[N];
    float * y = new float[N];
    float * y_reference = new float[N];

    std::cout << std::boolalpha;

    // Initialize x
    for (std::size_t i {0}; i < N; ++i) {
        x[i] = 0.1* (i % 10);
    }

    // Compute sin using the math library, to compare with my results
    for (size_t i {0}; i < N; ++i) {
        y_reference[i] = std::sin(x[i]);
    }

    // Variables to measure time
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    // Compute sin(x[i]) for each i, and output the results to y
    begin = std::chrono::steady_clock::now();
    lec2::sinx(N, terms, x, y);
    end = std::chrono::steady_clock::now();
    
    std::cout << "time sinx: " << std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count() << std::endl;
    std::cout << "Close to reference? " << lec2::are_close(N, y, y_reference, tolerance) << std::endl;

    // Compute again, but with 2 threads
    begin = std::chrono::steady_clock::now();
    lec2::parallel_sinx(N, terms, x, y);
    end = std::chrono::steady_clock::now();

    std::cout << "\ntime parallel_sinx: " << std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count() << std::endl;
    std::cout << "Close to reference? " << lec2::are_close(N, y, y_reference, tolerance) << std::endl;

    //***************************************** SIMD ******************************************
    //********************************* BAD ******************************************************
    std::cout << "****************** Bad init **********************" << std::endl;
    lec2::simd_example();

    lec2::bad_init_for_simd(N, x);

    begin = std::chrono::steady_clock::now();
    lec2::no_simd_baseline(N, x, y_reference);
    end = std::chrono::steady_clock::now();

    std::cout << "\ntime no SIMD: " << std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count() << std::endl;
    //print_array(8, y_reference);

    begin = std::chrono::steady_clock::now();
    lec2::simd_with_branching_AVX(N, x, y);
    end = std::chrono::steady_clock::now();

    std::cout << "\ntime SIMD: " << std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count() << std::endl;
    std::cout << "Equal reference? " << lec2::are_close(N, y, y_reference, 0) << std::endl;
    //print_array(8, y);
    
    //********************************* GOOD ******************************************************
    std::cout << "****************** Good init **********************" << std::endl;

    lec2::good_init_for_simd(N, x);

    begin = std::chrono::steady_clock::now();
    lec2::no_simd_baseline(N, x, y_reference);
    end = std::chrono::steady_clock::now();

    auto sequential_time {std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count()};
    std::cout << "\ntime no SIMD: " << sequential_time << std::endl;
    //print_array(8, y_reference);

    begin = std::chrono::steady_clock::now();
    lec2::simd_with_branching_AVX(N, x, y);
    end = std::chrono::steady_clock::now();

    auto not_optimized_simd_time {std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count()};
    std::cout << "\ntime SIMD: " << not_optimized_simd_time << std::endl;
    std::cout << "Equal reference? " << lec2::are_close(N, y, y_reference, 0) << std::endl;
    //print_array(8, y);

    begin = std::chrono::steady_clock::now();
    lec2::simd_optimized_for_good_init(N, x, y);
    end = std::chrono::steady_clock::now();

    auto optimized_simd_time {std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count()};
    std::cout << "\ntime SIMD optimized: " << optimized_simd_time << std::endl;
    std::cout << "Equal reference? " << lec2::are_close(N, y, y_reference, 0) << std::endl;
    //print_array(8, y);

    std::cout << "(sequential time)/(optimized simd time) = " << ((float)sequential_time)/optimized_simd_time << std::endl;
    std::cout << "(not optimized simd time)/(optimized simd time) = " << ((float)not_optimized_simd_time)/optimized_simd_time << std::endl;

    // Cleanup
    delete[] y_reference;
    delete[] y;
    delete[] x;
    y_reference = nullptr;
    y = nullptr;
    x = nullptr;
}

}