#include "common.h"

void bind_cuda_3d(py::module& m, py::module& s) {
    auto vector3m_class = bind<Vector3mC>(m, s, "Vector3m");
    auto vector3i_class = bind<Vector3iC>(m, s, "Vector3i");
    auto vector3u_class = bind<Vector3uC>(m, s, "Vector3u");
    auto vector3f_class = bind<Vector3fC>(m, s, "Vector3f");
    auto vector3d_class = bind<Vector3dC>(m, s, "Vector3d");

    vector3f_class
        .def(py::init<const Vector3f  &>())
        .def(py::init<const Vector3dC &>())
        .def(py::init<const Vector3iC &>())
        .def(py::init<const Vector3uC &>());

    vector3d_class
        .def(py::init<const Vector3d  &>())
        .def(py::init<const Vector3fC &>())
        .def(py::init<const Vector3iC &>())
        .def(py::init<const Vector3uC &>());

    vector3i_class
        .def(py::init<const Vector3i  &>())
        .def(py::init<const Vector3uC &>())
        .def(py::init<const Vector3fC &>())
        .def(py::init<const Vector3dC &>());

    vector3u_class
        .def(py::init<const Vector3u  &>())
        .def(py::init<const Vector3iC &>())
        .def(py::init<const Vector3fC &>())
        .def(py::init<const Vector3dC &>());
}
