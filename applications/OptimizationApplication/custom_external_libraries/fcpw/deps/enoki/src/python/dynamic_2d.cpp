#include "common.h"

void bind_dynamic_2d(py::module& m, py::module& s) {
    auto vector2m_class = bind<Vector2mX>(m, s, "Vector2m");
    auto vector2i_class = bind<Vector2iX>(m, s, "Vector2i");
    auto vector2u_class = bind<Vector2uX>(m, s, "Vector2u");
    auto vector2f_class = bind<Vector2fX>(m, s, "Vector2f");
    auto vector2d_class = bind<Vector2dX>(m, s, "Vector2d");

    vector2f_class
        .def(py::init<const Vector2f  &>())
        .def(py::init<const Vector2dX &>())
        .def(py::init<const Vector2uX &>())
        .def(py::init<const Vector2iX &>());

    vector2d_class
        .def(py::init<const Vector2d  &>())
        .def(py::init<const Vector2fX &>())
        .def(py::init<const Vector2uX &>())
        .def(py::init<const Vector2iX &>());

    vector2i_class
        .def(py::init<const Vector2i  &>())
        .def(py::init<const Vector2fX &>())
        .def(py::init<const Vector2dX &>())
        .def(py::init<const Vector2uX &>());

    vector2u_class
        .def(py::init<const Vector2u  &>())
        .def(py::init<const Vector2fX &>())
        .def(py::init<const Vector2dX &>())
        .def(py::init<const Vector2iX &>());
}
