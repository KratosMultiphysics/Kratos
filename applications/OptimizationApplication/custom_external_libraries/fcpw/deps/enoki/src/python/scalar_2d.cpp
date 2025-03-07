#include "common.h"

void bind_scalar_2d(py::module& m, py::module& s) {
    auto vector2m_class = bind<Vector2m>(m, s, "Vector2m");
    auto vector2i_class = bind<Vector2i>(m, s, "Vector2i");
    auto vector2u_class = bind<Vector2u>(m, s, "Vector2u");
    auto vector2f_class = bind<Vector2f>(m, s, "Vector2f");
    auto vector2d_class = bind<Vector2d>(m, s, "Vector2d");

    vector2f_class
        .def(py::init<const Vector2d &>())
        .def(py::init<const Vector2u &>())
        .def(py::init<const Vector2i &>());

    vector2d_class
        .def(py::init<const Vector2f &>())
        .def(py::init<const Vector2u &>())
        .def(py::init<const Vector2i &>());

    vector2i_class
        .def(py::init<const Vector2f &>())
        .def(py::init<const Vector2d &>())
        .def(py::init<const Vector2u &>());

    vector2u_class
        .def(py::init<const Vector2f &>())
        .def(py::init<const Vector2d &>())
        .def(py::init<const Vector2i &>());
}
