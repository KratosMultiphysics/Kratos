#include "common.h"
#include <pybind11/functional.h>

void bind_scalar_1d(py::module& m, py::module& s) {
    auto vector1m_class = bind<Vector1m>(m, s, "Vector1m");
    auto vector1i_class = bind<Vector1i>(m, s, "Vector1i");
    auto vector1u_class = bind<Vector1u>(m, s, "Vector1u");
    auto vector1f_class = bind<Vector1f>(m, s, "Vector1f");
    auto vector1d_class = bind<Vector1d>(m, s, "Vector1d");

    vector1f_class
        .def(py::init<const Vector1d &>())
        .def(py::init<const Vector1u &>())
        .def(py::init<const Vector1i &>());

    vector1d_class
        .def(py::init<const Vector1f &>())
        .def(py::init<const Vector1u &>())
        .def(py::init<const Vector1i &>());

    vector1i_class
        .def(py::init<const Vector1f &>())
        .def(py::init<const Vector1d &>())
        .def(py::init<const Vector1u &>());

    vector1u_class
        .def(py::init<const Vector1f &>())
        .def(py::init<const Vector1d &>())
        .def(py::init<const Vector1i &>());

    m.def(
        "binary_search",
        [](uint32_t start,
           uint32_t end,
           const std::function<bool(uint32_t)> &pred) {
            return enoki::binary_search(start, end, pred);
        },
        "start"_a, "end"_a, "pred"_a);
}
