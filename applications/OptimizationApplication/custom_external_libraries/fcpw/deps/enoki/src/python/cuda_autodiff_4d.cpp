#include "common.h"

void bind_cuda_autodiff_4d(py::module& m, py::module& s) {
    auto vector4m_class = bind<Vector4mD>(m, s, "Vector4m");
    auto vector4i_class = bind<Vector4iD>(m, s, "Vector4i");
    auto vector4u_class = bind<Vector4uD>(m, s, "Vector4u");
    auto vector4f_class = bind<Vector4fD>(m, s, "Vector4f");
    auto vector4d_class = bind<Vector4dD>(m, s, "Vector4d");

    vector4f_class
        .def(py::init<const Vector4f  &>())
        .def(py::init<const Vector4fC &>())
        .def(py::init<const Vector4dD &>())
        .def(py::init<const Vector4uD &>())
        .def(py::init<const Vector4iD &>());

    vector4d_class
        .def(py::init<const Vector4d  &>())
        .def(py::init<const Vector4dC &>())
        .def(py::init<const Vector4fD &>())
        .def(py::init<const Vector4uD &>())
        .def(py::init<const Vector4iD &>());

    vector4i_class
        .def(py::init<const Vector4i  &>())
        .def(py::init<const Vector4iC &>())
        .def(py::init<const Vector4fD &>())
        .def(py::init<const Vector4dD &>())
        .def(py::init<const Vector4uD &>());

    vector4u_class
        .def(py::init<const Vector4u  &>())
        .def(py::init<const Vector4uC &>())
        .def(py::init<const Vector4fD &>())
        .def(py::init<const Vector4dD &>())
        .def(py::init<const Vector4iD &>());
}
