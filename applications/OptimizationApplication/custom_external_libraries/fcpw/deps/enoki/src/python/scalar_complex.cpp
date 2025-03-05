#include "complex.h"

void bind_scalar_complex(py::module& m, py::module& s) {
    bind_complex<Complex2f>(m, s, "Complex2f");
    bind_complex<Complex24f>(m, s, "Complex24f");

    bind_complex<Complex2d>(m, s, "Complex2d");
    bind_complex<Complex24d>(m, s, "Complex24d");
}
