#include "common.h"

void bind_cuda_autodiff_0d(py::module& m, py::module& s) {
    auto vector0m_class = bind<Vector0mD>(m, s, "Vector0m");
    auto vector0i_class = bind<Vector0iD>(m, s, "Vector0i");
    auto vector0u_class = bind<Vector0uD>(m, s, "Vector0u");
    auto vector0f_class = bind<Vector0fD>(m, s, "Vector0f");
    auto vector0d_class = bind<Vector0dD>(m, s, "Vector0d");

    vector0f_class
        .def(py::init<const Vector0f  &>())
        .def(py::init<const Vector0fC &>())
        .def(py::init<const Vector0dD &>())
        .def(py::init<const Vector0uD &>())
        .def(py::init<const Vector0iD &>());

    vector0d_class
        .def(py::init<const Vector0d  &>())
        .def(py::init<const Vector0dC &>())
        .def(py::init<const Vector0fD &>())
        .def(py::init<const Vector0uD &>())
        .def(py::init<const Vector0iD &>());

    vector0i_class
        .def(py::init<const Vector0i  &>())
        .def(py::init<const Vector0iC &>())
        .def(py::init<const Vector0fD &>())
        .def(py::init<const Vector0dD &>())
        .def(py::init<const Vector0uD &>());

    vector0u_class
        .def(py::init<const Vector0u  &>())
        .def(py::init<const Vector0uC &>())
        .def(py::init<const Vector0fD &>())
        .def(py::init<const Vector0dD &>())
        .def(py::init<const Vector0iD &>());
}
