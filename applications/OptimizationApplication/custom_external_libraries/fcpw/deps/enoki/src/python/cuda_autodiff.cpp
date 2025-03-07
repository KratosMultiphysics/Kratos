#include "common.h"

extern void bind_cuda_autodiff_0d(py::module&, py::module&);
extern void bind_cuda_autodiff_1d(py::module&, py::module&);
extern void bind_cuda_autodiff_2d(py::module&, py::module&);
extern void bind_cuda_autodiff_3d(py::module&, py::module&);
extern void bind_cuda_autodiff_4d(py::module&, py::module&);
extern void bind_cuda_autodiff_complex(py::module&, py::module&);
extern void bind_cuda_autodiff_matrix(py::module&, py::module&);

bool *implicit_conversion = nullptr;

PYBIND11_MODULE(cuda_autodiff, s) {
    py::module m = py::module::import("enoki");
    py::module::import("enoki.cuda");

    implicit_conversion = (bool *) py::get_shared_data("implicit_conversion");

    bind_cuda_autodiff_1d(m, s);
    bind_cuda_autodiff_0d(m, s); // after FloatD
    bind_cuda_autodiff_2d(m, s);
    bind_cuda_autodiff_3d(m, s);
    bind_cuda_autodiff_4d(m, s);
    bind_cuda_autodiff_complex(m, s);
    bind_cuda_autodiff_matrix(m, s);

    m.def("set_requires_gradient",
          [](py::object o, bool value) {
              throw py::type_error("set_requires_gradient(): requires a differentiable type as input!");
          }, "array"_a, "value"_a = true);
}
