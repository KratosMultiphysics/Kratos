#include "common.h"

extern void bind_dynamic_0d(py::module&, py::module&);
extern void bind_dynamic_1d(py::module&, py::module&);
extern void bind_dynamic_2d(py::module&, py::module&);
extern void bind_dynamic_3d(py::module&, py::module&);
extern void bind_dynamic_4d(py::module&, py::module&);
extern void bind_dynamic_complex(py::module&, py::module&);
extern void bind_dynamic_matrix(py::module&, py::module&);
extern void bind_dynamic_pcg32(py::module&, py::module&);

bool *implicit_conversion = nullptr;

PYBIND11_MODULE(dynamic, s) {
    py::module m = py::module::import("enoki");
    py::module::import("enoki.scalar");

    implicit_conversion = (bool *) py::get_shared_data("implicit_conversion");

    bind_dynamic_1d(m, s);
    bind_dynamic_0d(m, s); // after FloatX
    bind_dynamic_2d(m, s);
    bind_dynamic_3d(m, s);
    bind_dynamic_4d(m, s);
    bind_dynamic_complex(m, s);
    bind_dynamic_matrix(m, s);
    bind_dynamic_pcg32(m, s);
}
