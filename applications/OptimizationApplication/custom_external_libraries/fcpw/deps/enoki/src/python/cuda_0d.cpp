#include "common.h"

void bind_cuda_0d(py::module& m, py::module& s) {
    auto vector0m_class = bind<Vector0mC>(m, s, "Vector0m");
    auto vector0i_class = bind<Vector0iC>(m, s, "Vector0i");
    auto vector0u_class = bind<Vector0uC>(m, s, "Vector0u");
    auto vector0f_class = bind<Vector0fC>(m, s, "Vector0f");
    auto vector0d_class = bind<Vector0dC>(m, s, "Vector0d");

    vector0f_class
        .def(py::init<const Vector0f  &>())
        .def(py::init<const Vector0dC &>())
        .def(py::init<const Vector0uC &>())
        .def(py::init<const Vector0iC &>());

    vector0d_class
        .def(py::init<const Vector0d  &>())
        .def(py::init<const Vector0fC &>())
        .def(py::init<const Vector0uC &>())
        .def(py::init<const Vector0iC &>());

    vector0i_class
        .def(py::init<const Vector0i  &>())
        .def(py::init<const Vector0fC &>())
        .def(py::init<const Vector0dC &>())
        .def(py::init<const Vector0uC &>());

    vector0u_class
        .def(py::init<const Vector0u  &>())
        .def(py::init<const Vector0fC &>())
        .def(py::init<const Vector0dC &>())
        .def(py::init<const Vector0iC &>());
}
