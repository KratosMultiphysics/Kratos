#include "common.h"

void bind_scalar_3d(py::module& m, py::module& s) {
    auto vector3m_class = bind<Vector3m>(m, s, "Vector3m");
    auto vector3i_class = bind<Vector3i>(m, s, "Vector3i");
    auto vector3u_class = bind<Vector3u>(m, s, "Vector3u");
    auto vector3f_class = bind<Vector3f>(m, s, "Vector3f");
    auto vector3d_class = bind<Vector3d>(m, s, "Vector3d");

    vector3f_class
        .def(py::init<const Vector3d &>())
        .def(py::init<const Vector3u &>())
        .def(py::init<const Vector3i &>());

    vector3d_class
        .def(py::init<const Vector3f &>())
        .def(py::init<const Vector3u &>())
        .def(py::init<const Vector3i &>());

    vector3i_class
        .def(py::init<const Vector3f &>())
        .def(py::init<const Vector3d &>())
        .def(py::init<const Vector3u &>());

    vector3u_class
        .def(py::init<const Vector3f &>())
        .def(py::init<const Vector3d &>())
        .def(py::init<const Vector3i &>());
}
