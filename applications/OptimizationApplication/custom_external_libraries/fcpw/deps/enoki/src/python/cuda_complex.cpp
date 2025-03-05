#include "complex.h"

void bind_cuda_complex(py::module& m, py::module& s) {
    bind_complex<Complex2fC>(m, s, "Complex2f")
        .def(py::init<const Complex2f &>());

    bind_complex<Complex24fC>(m, s, "Complex24f")
        .def(py::init<const Complex24f &>());

    bind_complex<Complex2dC>(m, s, "Complex2d")
        .def(py::init<const Complex2d &>());

    bind_complex<Complex24dC>(m, s, "Complex24d")
        .def(py::init<const Complex24d &>());
}
