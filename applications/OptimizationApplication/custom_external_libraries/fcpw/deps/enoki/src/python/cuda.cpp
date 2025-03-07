#include "common.h"

extern void bind_cuda_0d(py::module&, py::module&);
extern void bind_cuda_1d(py::module&, py::module&);
extern void bind_cuda_2d(py::module&, py::module&);
extern void bind_cuda_3d(py::module&, py::module&);
extern void bind_cuda_4d(py::module&, py::module&);
extern void bind_cuda_complex(py::module&, py::module&);
extern void bind_cuda_matrix(py::module&, py::module&);
extern void bind_cuda_pcg32(py::module&, py::module&);

bool *implicit_conversion = nullptr;

PYBIND11_MODULE(cuda, s) {
    py::module m = py::module::import("enoki");
    py::module::import("enoki.scalar");

    implicit_conversion = (bool *) py::get_shared_data("implicit_conversion");

    py::class_<Buffer<true>>(m, "GPUBuffer");

    cuda_sync();
    bind_cuda_1d(m, s);
    bind_cuda_0d(m, s); // after FloatC
    bind_cuda_2d(m, s);
    bind_cuda_3d(m, s);
    bind_cuda_4d(m, s);
    bind_cuda_complex(m, s);
    bind_cuda_matrix(m, s);
    bind_cuda_pcg32(m, s);

    m.def("cuda_eval", &cuda_eval, "log_assembly"_a = false,
          py::call_guard<py::gil_scoped_release>());

    m.def("cuda_sync", &cuda_sync,
          py::call_guard<py::gil_scoped_release>());

    m.def("cuda_malloc_trim", &cuda_malloc_trim);

    m.def("cuda_whos", []() { char *w = cuda_whos(); py::print(w); free(w); });

    m.def("cuda_mem_get_info", []() {
        size_t free = 0, total = 0;
        cuda_mem_get_info(&free, &total);
        return std::make_pair(free, total);
    });

    m.def("cuda_set_log_level", &cuda_set_log_level,
          "Sets the current log level (0: none, 1: kernel launches, 2: +ptxas "
          "statistics, 3: +ptx source, 4: +jit trace, 5: +ref counting)");

    m.def("cuda_log_level", &cuda_log_level);
}
