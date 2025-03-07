#include "common.h"

void bind_cuda_2d(py::module& m, py::module& s) {
    auto vector2m_class = bind<Vector2mC>(m, s, "Vector2m");
    auto vector2i_class = bind<Vector2iC>(m, s, "Vector2i");
    auto vector2u_class = bind<Vector2uC>(m, s, "Vector2u");
    auto vector2f_class = bind<Vector2fC>(m, s, "Vector2f");
    auto vector2d_class = bind<Vector2dC>(m, s, "Vector2d");

    vector2f_class
        .def(py::init<const Vector2f  &>())
        .def(py::init<const Vector2dC &>())
        .def(py::init<const Vector2uC &>())
        .def(py::init<const Vector2iC &>());

    vector2d_class
        .def(py::init<const Vector2d  &>())
        .def(py::init<const Vector2fC &>())
        .def(py::init<const Vector2uC &>())
        .def(py::init<const Vector2iC &>());

    vector2i_class
        .def(py::init<const Vector2i  &>())
        .def(py::init<const Vector2fC &>())
        .def(py::init<const Vector2dC &>())
        .def(py::init<const Vector2uC &>());

    vector2u_class
        .def(py::init<const Vector2u  &>())
        .def(py::init<const Vector2fC &>())
        .def(py::init<const Vector2dC &>())
        .def(py::init<const Vector2iC &>());
}
