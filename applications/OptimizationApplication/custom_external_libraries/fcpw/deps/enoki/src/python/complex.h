#pragma once
#include "common.h"
#include <enoki/complex.h>

template <typename Complex>
py::class_<Complex> bind_complex(py::module &m, py::module &s, const char *name) {
    using Value  = value_t<Complex>;
    using Scalar = scalar_t<Value>;
    using Mask   = mask_t<Complex>;

    auto cls = py::class_<Complex>(s, name)
        .def(py::init<>())
        .def(py::init<const Scalar &>());
        if constexpr (!std::is_same_v<Value, Scalar>)
            cls.def(py::init<const Value &>());
        cls.def(py::init<const Value &, const Value &>(), "real"_a, "imag"_a)
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def(py::self - py::self)
        .def(py::self + py::self)
        .def(py::self * py::self)
        .def(py::self / py::self)
        .def(-py::self)
        .def("__repr__", [](const Complex &a) -> std::string {
            if (*implicit_conversion)
                return "";
            std::ostringstream oss;
            oss << a;
            return oss.str();
        })
        .def("__getitem__", [](const Complex &a, size_t index) {
            if (index >= 2)
                throw py::index_error();
            return a.coeff(index);
        })
        .def("__setitem__", [](Complex &a, size_t index, const Value &value) {
            if (index >= 2)
                throw py::index_error();
            a.coeff(index) = value;
        })
        .def_static("identity", [](size_t size) { return identity<Complex>(size); }, "size"_a = 1)
        .def_static("zero", [](size_t size) { return zero<Complex>(size); }, "size"_a = 1)
        .def_static("full", [](Scalar value, size_t size) { return full<Complex>(value, size); },
                    "value"_a, "size"_a = 1);

    m.def("real", [](const Complex &z) { return real(z); });
    m.def("imag", [](const Complex &z) { return imag(z); });
    m.def("norm", [](const Complex &z) { return norm(z); });
    m.def("squared_norm", [](const Complex &z) { return squared_norm(z); });
    m.def("rcp", [](const Complex &z) { return rcp(z); });
    m.def("conj", [](const Complex &z) { return conj(z); });
    m.def("exp", [](const Complex &z) { return exp(z); });
    m.def("log", [](const Complex &z) { return log(z); });
    m.def("arg", [](const Complex &z) { return arg(z); });
    m.def("pow", [](const Complex &z1, const Complex &z2) { return pow(z1, z2); });
    m.def("sqrt", [](const Complex &z) { return sqrt(z); });
    m.def("sin", [](const Complex &z) { return sin(z); });
    m.def("cos", [](const Complex &z) { return cos(z); });
    m.def("sincos", [](const Complex &z) { return sincos(z); });
    m.def("tan", [](const Complex &z) { return tan(z); });
    m.def("asin", [](const Complex &z) { return asin(z); });
    m.def("acos", [](const Complex &z) { return acos(z); });
    m.def("atan", [](const Complex &z) { return atan(z); });
    m.def("sinh", [](const Complex &z) { return sinh(z); });
    m.def("cosh", [](const Complex &z) { return cosh(z); });
    m.def("sincosh", [](const Complex &z) { return sincosh(z); });
    m.def("tanh", [](const Complex &z) { return tanh(z); });
    m.def("asinh", [](const Complex &z) { return asinh(z); });
    m.def("acosh", [](const Complex &z) { return acosh(z); });
    m.def("atanh", [](const Complex &z) { return atanh(z); });

    m.def("isfinite", [](const Complex &z) -> Mask { return enoki::isfinite(z); });
    m.def("isnan",    [](const Complex &z) -> Mask { return enoki::isnan(z); });
    m.def("isinf",    [](const Complex &z) -> Mask { return enoki::isinf(z); });

    if constexpr (is_diff_array_v<Complex>) {
        using Detached = expr_t<decltype(detach(std::declval<Complex&>()))>;

        m.def("detach", [](const Complex &a) -> Detached { return detach(a); });
        m.def("requires_gradient",
              [](const Complex &a) { return requires_gradient(a); },
              "array"_a);

        m.def("set_requires_gradient",
              [](Complex &a, bool value) { set_requires_gradient(a, value); },
              "array"_a, "value"_a = true);

        m.def("gradient", [](Complex &a) { return eval(gradient(a)); });
        m.def("set_gradient",
              [](Complex &a, const Detached &g, bool b) { set_gradient(a, g, b); },
              "array"_a, "gradient"_a, "backward"_a = true);

        m.def("graphviz", [](const Complex &a) { return graphviz(a); });

        m.def("set_label", [](const Complex &a, const char *label) {
            set_label(a, label);
        });
    }

    implicitly_convertible<Value, Complex>();
    if constexpr (!std::is_same_v<Scalar, Value>)
        implicitly_convertible<Scalar, Complex>();

    return cls;
}
