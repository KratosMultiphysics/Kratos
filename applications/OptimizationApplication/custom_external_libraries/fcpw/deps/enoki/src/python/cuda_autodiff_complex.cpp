#include "complex.h"

void bind_cuda_autodiff_complex(py::module& m, py::module& s) {
    bind_complex<Complex2fD>(m, s, "Complex2f")
        .def(py::init<const Complex2f &>())
        .def(py::init<const Complex2fC &>());

    bind_complex<Complex24fD>(m, s, "Complex24f")
        .def(py::init<const Complex24f &>())
        .def(py::init<const Complex24fC &>());

    bind_complex<Complex2dD>(m, s, "Complex2d")
        .def(py::init<const Complex2d &>())
        .def(py::init<const Complex2dC &>());

    bind_complex<Complex24dD>(m, s, "Complex24d")
        .def(py::init<const Complex24d &>())
        .def(py::init<const Complex24dC &>());
}
