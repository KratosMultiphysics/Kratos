#include "common.h"

void bind_scalar_0d(py::module& m, py::module& s) {
    auto vector0m_class = bind<Vector0m>(m, s, "Vector0m");
    auto vector0i_class = bind<Vector0i>(m, s, "Vector0i");
    auto vector0u_class = bind<Vector0u>(m, s, "Vector0u");
    auto vector0f_class = bind<Vector0f>(m, s, "Vector0f");
    auto vector0d_class = bind<Vector0d>(m, s, "Vector0d");

    vector0f_class
        .def(py::init<const Vector0d &>())
        .def(py::init<const Vector0u &>())
        .def(py::init<const Vector0i &>());

    vector0d_class
        .def(py::init<const Vector0f &>())
        .def(py::init<const Vector0u &>())
        .def(py::init<const Vector0i &>());

    vector0i_class
        .def(py::init<const Vector0f &>())
        .def(py::init<const Vector0d &>())
        .def(py::init<const Vector0u &>());

    vector0u_class
        .def(py::init<const Vector0f &>())
        .def(py::init<const Vector0d &>())
        .def(py::init<const Vector0i &>());
}
