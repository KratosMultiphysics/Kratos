#include "common.h"
#include <pybind11/functional.h>

void bind_cuda_1d(py::module& m, py::module& s) {
    auto mask_class = bind<MaskC>(m, s, "Mask");
    auto uint32_class = bind<UInt32C>(m, s, "UInt32");
    auto uint64_class = bind<UInt64C>(m, s, "UInt64");
    auto int32_class = bind<Int32C>(m, s, "Int32");
    auto int64_class = bind<Int64C>(m, s, "Int64");
    auto float32_class = bind<Float32C>(m, s, "Float32");
    auto float64_class = bind<Float64C>(m, s, "Float64");

    float32_class
        .def(py::init<const Float64C &>())
        .def(py::init<const Int32C &>())
        .def(py::init<const Int64C &>())
        .def(py::init<const UInt32C &>())
        .def(py::init<const UInt64C &>());

    float64_class
        .def(py::init<const Float32C &>())
        .def(py::init<const Int32C &>())
        .def(py::init<const Int64C &>())
        .def(py::init<const UInt32C &>())
        .def(py::init<const UInt64C &>());

    int32_class
        .def(py::init<const Float32C &>())
        .def(py::init<const Float64C &>())
        .def(py::init<const Int64C &>())
        .def(py::init<const UInt32C &>())
        .def(py::init<const UInt64C &>());

    int64_class
        .def(py::init<const Float32C &>())
        .def(py::init<const Float64C &>())
        .def(py::init<const Int32C &>())
        .def(py::init<const UInt32C &>())
        .def(py::init<const UInt64C &>());

    uint32_class
        .def(py::init<const Float32C &>())
        .def(py::init<const Float64C &>())
        .def(py::init<const Int32C &>())
        .def(py::init<const Int64C &>())
        .def(py::init<const UInt64C &>());

    uint64_class
        .def(py::init<const Float32C &>())
        .def(py::init<const Float64C &>())
        .def(py::init<const Int32C &>())
        .def(py::init<const Int64C &>())
        .def(py::init<const UInt32C &>());

    auto vector1m_class = bind<Vector1mC>(m, s, "Vector1m");
    auto vector1i_class = bind<Vector1iC>(m, s, "Vector1i");
    auto vector1u_class = bind<Vector1uC>(m, s, "Vector1u");
    auto vector1f_class = bind<Vector1fC>(m, s, "Vector1f");
    auto vector1d_class = bind<Vector1dC>(m, s, "Vector1d");

    vector1f_class
        .def(py::init<const Vector1f  &>())
        .def(py::init<const Vector1dC &>())
        .def(py::init<const Vector1uC &>())
        .def(py::init<const Vector1iC &>());

    vector1d_class
        .def(py::init<const Vector1d  &>())
        .def(py::init<const Vector1fC &>())
        .def(py::init<const Vector1uC &>())
        .def(py::init<const Vector1iC &>());

    vector1i_class
        .def(py::init<const Vector1i  &>())
        .def(py::init<const Vector1fC &>())
        .def(py::init<const Vector1dC &>())
        .def(py::init<const Vector1uC &>());

    vector1u_class
        .def(py::init<const Vector1u  &>())
        .def(py::init<const Vector1fC &>())
        .def(py::init<const Vector1dC &>())
        .def(py::init<const Vector1iC &>());

    m.def(
        "binary_search",
        [](uint32_t start,
           uint32_t end,
           const std::function<MaskC(UInt32C)> &pred) {
            return enoki::binary_search(start, end, pred);
        },
        "start"_a, "end"_a, "pred"_a);

    m.def("meshgrid", [](const Float32C &x, const Float32C &y) {
        auto result = meshgrid(x, y);
        return std::make_pair(std::move(result.x()), std::move(result.y()));
    });

    m.def("meshgrid", [](const Float64C &x, const Float64C &y) {
        auto result = meshgrid(x, y);
        return std::make_pair(std::move(result.x()), std::move(result.y()));
    });

    m.def("partition", [](const UInt64C &x) {
        return partition(x);
    });
}
