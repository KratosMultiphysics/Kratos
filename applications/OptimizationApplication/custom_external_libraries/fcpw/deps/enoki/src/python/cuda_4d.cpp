#include "common.h"

void bind_cuda_4d(py::module& m, py::module& s) {
    auto vector4m_class = bind<Vector4mC>(m, s, "Vector4m");
    auto vector4i_class = bind<Vector4iC>(m, s, "Vector4i");
    auto vector4u_class = bind<Vector4uC>(m, s, "Vector4u");
    auto vector4f_class = bind<Vector4fC>(m, s, "Vector4f");
    auto vector4d_class = bind<Vector4dC>(m, s, "Vector4d");

    vector4f_class
        .def(py::init<const Vector4f  &>())
        .def(py::init<const Vector4dC &>())
        .def(py::init<const Vector4iC &>())
        .def(py::init<const Vector4uC &>());

    vector4d_class
        .def(py::init<const Vector4d  &>())
        .def(py::init<const Vector4fC &>())
        .def(py::init<const Vector4iC &>())
        .def(py::init<const Vector4uC &>());

    vector4i_class
        .def(py::init<const Vector4i  &>())
        .def(py::init<const Vector4uC &>())
        .def(py::init<const Vector4fC &>())
        .def(py::init<const Vector4dC &>());

    vector4u_class
        .def(py::init<const Vector4u  &>())
        .def(py::init<const Vector4iC &>())
        .def(py::init<const Vector4fC &>())
        .def(py::init<const Vector4dC &>());
}
