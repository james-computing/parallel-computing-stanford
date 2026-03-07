#include "../../include/hpp/mandelbrot.hpp"

void run_mandelbrot_example() {
    unsigned int const width = 768;
    unsigned int const height= 512;

    float const x0 = -2.;
    float const x1 = 1.;
    float const y0 = -1.;
    float const y1 = 1.;

    int const maxIterations = 256;
    int * const buf = new int[width * height]; // output

    ispc::mandelbrot_ispc(x0, y0, x1, y1, width, height, maxIterations, buf);
}