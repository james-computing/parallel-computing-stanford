export void ispc_sinx_interleaved(uniform size_t const N, uniform size_t const terms, uniform float const x[], uniform float result[]) {
    // assumes N % programCount = 0
    for (uniform int i = 0; i < N; i += programCount) {
        int idx = i + programIndex;
        float value = x[idx];
        float number = x[idx] * x[idx] * x[idx];
        uniform int denom = 6; // 3!
        uniform int sign = -1;

        for (uniform int j = 1; j <= terms; ++j) {
            value += sign * number/denom;
            number *= x[idx] * x[idx];
            denom *= (2*j + 2)*(2*j + 3);
            sign *= -1;
        }
        result[idx] = value;
    }
}

export void ispc_sinx_blocks(uniform size_t const N, uniform size_t const terms, uniform float const x[], uniform float result[]) {
    // assumes N % programCount = 0
    uniform int blockSize = N / programCount;
    int start = blockSize * programIndex;
    for (uniform int i = 0; i < blockSize; ++i) {
        int idx = start + i;
        float value = x[idx];
        float number = x[idx] * x[idx] * x[idx];
        uniform int denom = 6; // 3!
        uniform int sign = -1;

        for (uniform int j = 1; j <= terms; ++j) {
            value += sign * number/denom;
            number *= x[idx] * x[idx];
            denom *= (2*j + 2)*(2*j + 3);
            sign *= -1;
        }
        result[idx] = value;
    }
}

export void ispc_sinx_abstract(uniform size_t const N, uniform size_t const terms, uniform float const x[], uniform float result[]) {
    foreach (i = 0 ... N) {
        float value = x[i];
        float number = x[i] * x[i] * x[i];
        uniform int denom = 6; // 3!
        uniform int sign = -1;

        for (uniform int j = 1; j <= terms; ++j) {
            value += sign * number/denom;
            number *= x[i] * x[i];
            denom *= (2*j + 2)*(2*j + 3);
            sign *= -1;
        }
        result[i] = value;
    }
}