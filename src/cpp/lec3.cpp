#include "../../include/hpp/lec2.hpp"
#include "../../include/hpp/lec3.hpp"

namespace lec3 {

void execute_examples() {
    size_t const exponent {25};
    size_t const N {1 << exponent};
    size_t const terms {10};
    float * x {new float[N]};
    float * y {new float[N]};
    float * y_ref {new float[N]};

    float const tolerance {1e-7};
    bool closeToReference;

    // Variables to measure time
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    std::cout << std::boolalpha;

    lec2::init_for_sinx(N, x);
    
    for (size_t i {0}; i < N; ++i) {
        y_ref[i] = std::sin(x[i]);
    }

    std::cout << "\nsinx" << std::endl;
    begin = std::chrono::steady_clock::now();
    lec2::sinx(N, terms, x, y);
    end = std::chrono::steady_clock::now();

    closeToReference = lec2::are_close(N, y, y_ref, tolerance);
    std::cout << "Close to reference? " << closeToReference << std::endl;
    std::cout << "Time = " << std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count() << std::endl;

    std::cout << "\nparallel_sinx" << std::endl;
    begin = std::chrono::steady_clock::now();
    lec2::parallel_sinx(N, terms, x, y);
    end = std::chrono::steady_clock::now();

    closeToReference = lec2::are_close(N, y, y_ref, tolerance);
    std::cout << "Equal to reference? " << closeToReference << std::endl;
    std::cout << "Time = " << std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count() << std::endl;

    std::cout << "\nispc_sinx_interleaved" << std::endl;
    begin = std::chrono::steady_clock::now();
    ispc::ispc_sinx_interleaved(N, terms, x, y);
    end = std::chrono::steady_clock::now();

    closeToReference = lec2::are_close(N, y, y_ref, tolerance);
    std::cout << "Equal to reference? " << closeToReference << std::endl;
    std::cout << "Time = " << std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count() << std::endl;

    std::cout << "\nispc_sinx_blocks" << std::endl;
    begin = std::chrono::steady_clock::now();
    ispc::ispc_sinx_blocks(N, terms, x, y);
    end = std::chrono::steady_clock::now();

    closeToReference = lec2::are_close(N, y, y_ref, tolerance);
    std::cout << "Equal to reference? " << closeToReference << std::endl;
    std::cout << "Time = " << std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count() << std::endl;

    std::cout << "\nispc_sinx_abstract" << std::endl;
    begin = std::chrono::steady_clock::now();
    ispc::ispc_sinx_abstract(N, terms, x, y);
    end = std::chrono::steady_clock::now();

    closeToReference = lec2::are_close(N, y, y_ref, tolerance);
    std::cout << "Equal to reference? " << closeToReference << std::endl;
    std::cout << "Time = " << std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count() << std::endl;

    // cleanup
    delete[] y_ref;
    delete[] y;
    delete[] x;
    y_ref = nullptr;
    y = nullptr;
    x = nullptr;
}

}