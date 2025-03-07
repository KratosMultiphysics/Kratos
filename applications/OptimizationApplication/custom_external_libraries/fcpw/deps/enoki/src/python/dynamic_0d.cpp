#include "common.h"

void bind_dynamic_0d(py::module& m, py::module& s) {
    auto vector0m_class = bind<Vector0mX>(m, s, "Vector0m");
    auto vector0i_class = bind<Vector0iX>(m, s, "Vector0i");
    auto vector0u_class = bind<Vector0uX>(m, s, "Vector0u");
    auto vector0f_class = bind<Vector0fX>(m, s, "Vector0f");
    auto vector0d_class = bind<Vector0dX>(m, s, "Vector0d");

    vector0f_class
        .def(py::init<const Vector0f  &>())
        .def(py::init<const Vector0dX &>())
        .def(py::init<const Vector0uX &>())
        .def(py::init<const Vector0iX &>());

    vector0d_class
        .def(py::init<const Vector0d  &>())
        .def(py::init<const Vector0fX &>())
        .def(py::init<const Vector0uX &>())
        .def(py::init<const Vector0iX &>());

    vector0i_class
        .def(py::init<const Vector0i  &>())
        .def(py::init<const Vector0fX &>())
        .def(py::init<const Vector0dX &>())
        .def(py::init<const Vector0uX &>());

    vector0u_class
        .def(py::init<const Vector0u  &>())
        .def(py::init<const Vector0fX &>())
        .def(py::init<const Vector0dX &>())
        .def(py::init<const Vector0iX &>());
}
