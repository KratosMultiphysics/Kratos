#include "common.h"

void bind_dynamic_3d(py::module& m, py::module& s) {
    auto vector3m_class = bind<Vector3mX>(m, s, "Vector3m");
    auto vector3i_class = bind<Vector3iX>(m, s, "Vector3i");
    auto vector3u_class = bind<Vector3uX>(m, s, "Vector3u");
    auto vector3f_class = bind<Vector3fX>(m, s, "Vector3f");
    auto vector3d_class = bind<Vector3dX>(m, s, "Vector3d");

    vector3f_class
        .def(py::init<const Vector3f  &>())
        .def(py::init<const Vector3dX &>())
        .def(py::init<const Vector3uX &>())
        .def(py::init<const Vector3iX &>());

    vector3d_class
        .def(py::init<const Vector3d  &>())
        .def(py::init<const Vector3fX &>())
        .def(py::init<const Vector3uX &>())
        .def(py::init<const Vector3iX &>());

    vector3i_class
        .def(py::init<const Vector3i  &>())
        .def(py::init<const Vector3fX &>())
        .def(py::init<const Vector3dX &>())
        .def(py::init<const Vector3uX &>());

    vector3u_class
        .def(py::init<const Vector3u  &>())
        .def(py::init<const Vector3fX &>())
        .def(py::init<const Vector3dX &>())
        .def(py::init<const Vector3iX &>());
}
