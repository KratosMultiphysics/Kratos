#include "random.h"

void bind_cuda_pcg32(py::module& m, py::module& s) {
    bind_pcg32<PCG32<Float32C, 1>>(m, s, "PCG32");
}
