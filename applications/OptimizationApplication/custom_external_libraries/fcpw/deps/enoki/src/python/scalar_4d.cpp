#include "common.h"


void bind_scalar_4d(py::module& m, py::module& s) {
    auto vector4m_class = bind<Vector4m>(m, s, "Vector4m");
    auto vector4i_class = bind<Vector4i>(m, s, "Vector4i");
    auto vector4u_class = bind<Vector4u>(m, s, "Vector4u");
    auto vector4f_class = bind<Vector4f>(m, s, "Vector4f");
    auto vector4d_class = bind<Vector4d>(m, s, "Vector4d");

    vector4f_class
        .def(py::init<const Vector4d &>())
        .def(py::init<const Vector4u &>())
        .def(py::init<const Vector4i &>());

    vector4d_class
        .def(py::init<const Vector4f &>())
        .def(py::init<const Vector4u &>())
        .def(py::init<const Vector4i &>());

    vector4i_class
        .def(py::init<const Vector4f &>())
        .def(py::init<const Vector4d &>())
        .def(py::init<const Vector4u &>());

    vector4u_class
        .def(py::init<const Vector4f &>())
        .def(py::init<const Vector4d &>())
        .def(py::init<const Vector4i &>());
}
