#include "random.h"

void bind_scalar_pcg32(py::module& m, py::module& s) {
    bind_pcg32<PCG32<float, 1>>(m, s, "PCG32");
}
