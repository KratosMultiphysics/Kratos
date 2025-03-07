#include "common.h"
#include <pybind11/functional.h>

void bind_dynamic_1d(py::module& m, py::module& s) {
    auto mask_class = bind<MaskX>(m, s, "Mask");
    auto mask64_class = bind<Mask64X>(m, s, "Mask64");
    auto uint32_class = bind<UInt32X>(m, s, "UInt32");
    auto uint64_class = bind<UInt64X>(m, s, "UInt64");
    auto int32_class = bind<Int32X>(m, s, "Int32");
    auto int64_class = bind<Int64X>(m, s, "Int64");
    auto float32_class = bind<Float32X>(m, s, "Float32");
    auto float64_class = bind<Float64X>(m, s, "Float64");

    mask_class
        .def(py::init<const Mask64X &>());

    mask64_class
        .def(py::init<const MaskX &>());

    implicitly_convertible<Mask64X, MaskX>();
    implicitly_convertible<MaskX, Mask64X>();

    float32_class
        .def(py::init<const Float64X &>())
        .def(py::init<const Int32X &>())
        .def(py::init<const Int64X &>())
        .def(py::init<const UInt32X &>())
        .def(py::init<const UInt64X &>());

    float64_class
        .def(py::init<const Float32X &>())
        .def(py::init<const Int32X &>())
        .def(py::init<const Int64X &>())
        .def(py::init<const UInt32X &>())
        .def(py::init<const UInt64X &>());

    int32_class
        .def(py::init<const Float32X &>())
        .def(py::init<const Float64X &>())
        .def(py::init<const Int64X &>())
        .def(py::init<const UInt32X &>())
        .def(py::init<const UInt64X &>());

    int64_class
        .def(py::init<const Float32X &>())
        .def(py::init<const Float64X &>())
        .def(py::init<const Int32X &>())
        .def(py::init<const UInt32X &>())
        .def(py::init<const UInt64X &>());

    uint32_class
        .def(py::init<const Float32X &>())
        .def(py::init<const Float64X &>())
        .def(py::init<const Int32X &>())
        .def(py::init<const Int64X &>())
        .def(py::init<const UInt64X &>());

    uint64_class
        .def(py::init<const Float32X &>())
        .def(py::init<const Float64X &>())
        .def(py::init<const Int32X &>())
        .def(py::init<const Int64X &>())
        .def(py::init<const UInt32X &>());

    auto vector1m_class = bind<Vector1mX>(m, s, "Vector1m");
    auto vector1i_class = bind<Vector1iX>(m, s, "Vector1i");
    auto vector1u_class = bind<Vector1uX>(m, s, "Vector1u");
    auto vector1f_class = bind<Vector1fX>(m, s, "Vector1f");
    auto vector1d_class = bind<Vector1dX>(m, s, "Vector1d");

    vector1f_class
        .def(py::init<const Vector1f  &>())
        .def(py::init<const Vector1dX &>())
        .def(py::init<const Vector1uX &>())
        .def(py::init<const Vector1iX &>());

    vector1d_class
        .def(py::init<const Vector1d  &>())
        .def(py::init<const Vector1fX &>())
        .def(py::init<const Vector1uX &>())
        .def(py::init<const Vector1iX &>());

    vector1i_class
        .def(py::init<const Vector1i  &>())
        .def(py::init<const Vector1fX &>())
        .def(py::init<const Vector1dX &>())
        .def(py::init<const Vector1uX &>());

    vector1u_class
        .def(py::init<const Vector1u  &>())
        .def(py::init<const Vector1fX &>())
        .def(py::init<const Vector1dX &>())
        .def(py::init<const Vector1iX &>());

    m.def(
        "binary_search",
        [](uint32_t start,
           uint32_t end,
           const std::function<MaskX (const UInt32X &)> &pred) {
            return enoki::binary_search(start, end, pred);
        },
        "start"_a, "end"_a, "pred"_a);

    m.def("meshgrid", [](const Float32X &x, const Float32X &y) {
        auto result = meshgrid(x, y);
        return std::make_pair(std::move(result.x()), std::move(result.y()));
    });

    m.def("meshgrid", [](const Float64X &x, const Float64X &y) {
        auto result = meshgrid(x, y);
        return std::make_pair(std::move(result.x()), std::move(result.y()));
    });
}
