#include "random.h"

void bind_dynamic_pcg32(py::module& m, py::module& s) {
    bind_pcg32<PCG32<Float32X, 1>>(m, s, "PCG32");
}
