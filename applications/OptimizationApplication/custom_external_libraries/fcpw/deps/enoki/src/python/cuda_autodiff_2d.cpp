#include "common.h"

void bind_cuda_autodiff_2d(py::module& m, py::module& s) {
    auto vector2m_class = bind<Vector2mD>(m, s, "Vector2m");
    auto vector2i_class = bind<Vector2iD>(m, s, "Vector2i");
    auto vector2u_class = bind<Vector2uD>(m, s, "Vector2u");
    auto vector2f_class = bind<Vector2fD>(m, s, "Vector2f");
    auto vector2d_class = bind<Vector2dD>(m, s, "Vector2d");

    vector2f_class
        .def(py::init<const Vector2f  &>())
        .def(py::init<const Vector2fC &>())
        .def(py::init<const Vector2dD &>())
        .def(py::init<const Vector2uD &>())
        .def(py::init<const Vector2iD &>());

    vector2d_class
        .def(py::init<const Vector2d  &>())
        .def(py::init<const Vector2dC &>())
        .def(py::init<const Vector2fD &>())
        .def(py::init<const Vector2uD &>())
        .def(py::init<const Vector2iD &>());

    vector2i_class
        .def(py::init<const Vector2i  &>())
        .def(py::init<const Vector2iC &>())
        .def(py::init<const Vector2fD &>())
        .def(py::init<const Vector2dD &>())
        .def(py::init<const Vector2uD &>());

    vector2u_class
        .def(py::init<const Vector2u  &>())
        .def(py::init<const Vector2uC &>())
        .def(py::init<const Vector2fD &>())
        .def(py::init<const Vector2dD &>())
        .def(py::init<const Vector2iD &>());
}
