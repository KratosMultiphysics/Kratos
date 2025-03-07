#include "complex.h"

void bind_dynamic_complex(py::module& m, py::module& s) {
    bind_complex<Complex2fX>(m, s, "Complex2f")
        .def(py::init<const Complex2f &>());

    bind_complex<Complex24fX>(m, s, "Complex24f")
        .def(py::init<const Complex24f &>());

    bind_complex<Complex2dX>(m, s, "Complex2d")
        .def(py::init<const Complex2d &>());

    bind_complex<Complex24dX>(m, s, "Complex24d")
        .def(py::init<const Complex24d &>());
}
