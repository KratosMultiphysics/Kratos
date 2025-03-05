#include "common.h"

void bind_dynamic_4d(py::module& m, py::module& s) {
    auto vector4m_class = bind<Vector4mX>(m, s, "Vector4m");
    auto vector4i_class = bind<Vector4iX>(m, s, "Vector4i");
    auto vector4u_class = bind<Vector4uX>(m, s, "Vector4u");
    auto vector4f_class = bind<Vector4fX>(m, s, "Vector4f");
    auto vector4d_class = bind<Vector4dX>(m, s, "Vector4d");

    vector4f_class
        .def(py::init<const Vector4f  &>())
        .def(py::init<const Vector4dX &>())
        .def(py::init<const Vector4uX &>())
        .def(py::init<const Vector4iX &>());

    vector4d_class
        .def(py::init<const Vector4d  &>())
        .def(py::init<const Vector4fX &>())
        .def(py::init<const Vector4uX &>())
        .def(py::init<const Vector4iX &>());

    vector4i_class
        .def(py::init<const Vector4i  &>())
        .def(py::init<const Vector4fX &>())
        .def(py::init<const Vector4dX &>())
        .def(py::init<const Vector4uX &>());

    vector4u_class
        .def(py::init<const Vector4u  &>())
        .def(py::init<const Vector4fX &>())
        .def(py::init<const Vector4dX &>())
        .def(py::init<const Vector4iX &>());
}
