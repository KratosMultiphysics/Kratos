#pragma once
#include "common.h"
#include <enoki/quaternion.h>

template <typename Quat>
py::class_<Quat> bind_quaternion(py::module &m, py::module &s, const char *name) {
    using Value  = value_t<Quat>;
    using Scalar = scalar_t<Quat>;
    using Mask   = mask_t<Quat>;

    auto cls = py::class_<Quat>(s, name)
        .def(py::init<>())
        .def(py::init<const Value &>(), "w"_a)
        .def(py::init<const Value &, const Value &, const Value &, const Value &>(), "x"_a, "y"_a, "z"_a, "w"_a)
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def(py::self - py::self)
        .def(py::self + py::self)
        .def(py::self * Value())
        .def(py::self / Value())
        .def(py::self * py::self)
        .def(py::self / py::self)
        .def(-py::self)
        .def("__repr__", [](const Quat &a) -> std::string {
            if (*implicit_conversion)
                return "";
            std::ostringstream oss;
            oss << a;
            return oss.str();
        })
        .def("__getitem__", [](const Quat &a, size_t index) {
            if (index >= 4)
                throw py::index_error();
            return a.coeff(index);
        })
        .def("__setitem__", [](Quat &a, size_t index, const Value &value) {
            if (index >= 4)
                throw py::index_error();
            a.coeff(index) = value;
        })
        .def("__setitem__", [](Quat &a, const mask_t<Value> &m, const Quat &b) {
            a[m] = b;
        })
        .def_static("identity", [](size_t size) { return identity<Quat>(size); }, "size"_a = 1)
        .def_static("zero", [](size_t size) { return zero<Quat>(size); }, "size"_a = 1)
        .def_static("full", [](Scalar value, size_t size) { return full<Quat>(value, size); },
                    "value"_a, "size"_a = 1);

    m.def("real", [](const Quat &a) { return real(a); });
    m.def("imag", [](const Quat &a) { return imag(a); });
    m.def("norm", [](const Quat &a) { return norm(a); });
    m.def("squared_norm", [](const Quat &a) { return squared_norm(a); });
    m.def("rcp", [](const Quat &a) { return rcp(a); });
    m.def("normalize", [](const Quat &a) { return normalize(a); });
    m.def("dot", [](const Quat &a, const Quat &b) { return dot(a, b); });

    m.def("abs", [](const Quat &a) { return abs(a); });
    m.def("sqrt", [](const Quat &a) { return sqrt(a); });
    m.def("exp", [](const Quat &a) { return exp(a); });
    m.def("log", [](const Quat &a) { return log(a); });
    m.def("pow", [](const Quat &a, const Quat &b) { return pow(a, b); });

    cls.def_property("x", [](const Quat &a) { return a.x(); },
                          [](Quat &a, const Value &v) { a.x() = v; });
    cls.def_property("y", [](const Quat &a) { return a.y(); },
                          [](Quat &a, const Value &v) { a.y() = v; });
    cls.def_property("z", [](const Quat &a) { return a.z(); },
                          [](Quat &a, const Value &v) { a.z() = v; });
    cls.def_property("w", [](const Quat &a) { return a.w(); },
                          [](Quat &a, const Value &v) { a.w() = v; });

    m.def("isfinite", [](const Quat &a) -> Mask { return enoki::isfinite(a); });
    m.def("isnan",    [](const Quat &a) -> Mask { return enoki::isnan(a); });
    m.def("isinf",    [](const Quat &a) -> Mask { return enoki::isinf(a); });

    using Vector3f = Array<Value, 3>;
    using Matrix4f = Matrix<Value, 4>;

    m.def("slerp",
          [](const Quat &a, const Quat &b, const Value &t) {
              return slerp(a, b, t);
          },
          "a"_a, "b"_a, "t"_a);

    m.def("quat_to_euler", [](const Quat &q) { return quat_to_euler<Vector3f>(q); });
    m.def("quat_to_matrix", [](const Quat &q) { return quat_to_matrix<Matrix4f>(q); });
    m.def("matrix_to_quat", [](const Matrix4f &m) { return matrix_to_quat(m); });

    m.def("rotate", [](const Vector3f &axis, const Value &angle) {
        return rotate<Quat>(axis, angle);
    }, "axis"_a, "angle"_a);

    return cls;
}
