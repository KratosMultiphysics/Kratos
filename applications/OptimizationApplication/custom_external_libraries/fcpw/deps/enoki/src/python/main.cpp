#include "common.h"
#include <pybind11/functional.h>

bool __implicit_conversion = false;

bool allclose_py(const py::object &a, const py::object &b,
                 const py::float_ &rtol, const py::float_ &atol,
                 bool equal_nan) {
    const char *tp_name_a = a.ptr()->ob_type->tp_name,
               *tp_name_b = b.ptr()->ob_type->tp_name;

    ssize_t la = PyObject_Length(a.ptr()),
            lb = PyObject_Length(b.ptr());

    bool num_a     = PyNumber_Check(a.ptr()) && la == -1,
         num_b     = PyNumber_Check(b.ptr()) && lb == -1,
         enoki_a   = strncmp(tp_name_a, "enoki.", 6) == 0,
         enoki_b   = strncmp(tp_name_b, "enoki.", 6) == 0,
         ndarray_a = strcmp(tp_name_a, "numpy.ndarray") == 0,
         ndarray_b = strcmp(tp_name_b, "numpy.ndarray") == 0;

    if (la == -1 || lb == -1)
        PyErr_Clear();

    if (enoki_a && (ndarray_b || num_b))
        return allclose_py(a, a.get_type()(b), rtol, atol, equal_nan);
    else if (enoki_b && (ndarray_a || num_a))
        return allclose_py(b.get_type()(a), b, rtol, atol, equal_nan);

    if (la != lb && !(((num_a || la == 1) && lb > 0) || ((num_b || lb == 1) && la > 0)))
        throw std::runtime_error("enoki.allclose(): length mismatch!");

    if ((enoki_a && enoki_b) || (num_a && num_b)) {
        py::module ek = py::module::import("enoki");

        py::object abs        = ek.attr("abs"),
                   eq         = ek.attr("eq"),
                   isnan      = ek.attr("isnan"),
                   isinf      = ek.attr("isinf"),
                   full       = ek.attr("full"),
                   all_nested = ek.attr("all_nested");

        py::object lhs = abs(a - b),
                   rhs = (num_b ? atol + abs(b) * rtol
                                : full(b.get_type(), atol) + abs(b) * rtol);

        py::object cond =
            py::reinterpret_steal<py::object>(PyObject_RichCompare(lhs.ptr(), rhs.ptr(), Py_LE));

        if (!cond)
            throw py::error_already_set();

        cond = cond | (isinf(a) & isinf(b));

        if (equal_nan)
            cond = cond | (isnan(a) & isnan(b));

        return py::cast<bool>(all_nested(cond));
    } else if (la >= 0) {
        for (size_t i = 0; i < (size_t) la; ++i) {
            py::int_ key(i);
            py::object ai = num_a ? a : a[key],
                       bi = num_b ? b : b[key];
            if (!allclose_py(ai, bi, rtol, atol, equal_nan))
                return false;
        }
    } else {
        throw std::runtime_error("enoki.allclose(): unsupported type!");
    }

    return true;
}

bool is_enoki_type(py::handle h) {
    return PyType_Check(h.ptr()) &&
           strncmp(((PyTypeObject *) h.ptr())->tp_name, "enoki.", 6) == 0;
}

PYBIND11_MODULE(core, m_) {
    ENOKI_MARK_USED(m_);
    py::module m = py::module::import("enoki");

    m.attr("__version__") = ENOKI_VERSION;
    py::set_shared_data("implicit_conversion", &__implicit_conversion);

    py::class_<Buffer<false>>(m, "CPUBuffer");

    m.def("empty",
        [](py::handle h, size_t size) {
            if (!is_enoki_type(h) && size == 1)
                return h();
            else
                return h.attr("empty")(size);
        },
        "type"_a, "size"_a = 1);

    m.def("zero",
        [](py::handle h, size_t size) {
            if (!is_enoki_type(h) && size == 1)
                return h(0);
            else
                return h.attr("zero")(size);
        },
        "type"_a, "size"_a = 1);

    m.def("arange",
        [](py::handle h, size_t size) {
            if (!is_enoki_type(h) && size == 1)
                return h(0);
            else
                return h.attr("arange")(size);
        },
        "type"_a, "size"_a = 1);

    m.def("full",
        [](py::handle h, py::handle value, size_t size) {
            if (!is_enoki_type(h) && size == 1)
                return h(value);
            else
                return h.attr("full")(value, size);
        },
        "type"_a, "value"_a, "size"_a = 1);

    m.def("linspace",
        [](py::handle h, py::handle start, py::handle end, size_t size) {
            if (!is_enoki_type(h))
                return h(start);
            else
                return h.attr("linspace")(start, end, size);
        },
        "type"_a, "start"_a, "end"_a, "size"_a = 1);

    m.def("allclose", &allclose_py,
        "a"_a, "b"_a, "rtol"_a = 1e-5, "atol"_a = 1e-8,
        "equal_nan"_a = false
    );

    m.attr("pi") = M_PI;
    m.attr("e") = M_E;
    m.attr("inf") = std::numeric_limits<float>::infinity();
    m.attr("nan") = std::numeric_limits<float>::quiet_NaN();
}
