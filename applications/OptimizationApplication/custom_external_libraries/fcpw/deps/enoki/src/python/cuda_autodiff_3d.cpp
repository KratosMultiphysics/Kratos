#include "common.h"

void bind_cuda_autodiff_3d(py::module& m, py::module& s) {
    auto vector3m_class = bind<Vector3mD>(m, s, "Vector3m");
    auto vector3i_class = bind<Vector3iD>(m, s, "Vector3i");
    auto vector3u_class = bind<Vector3uD>(m, s, "Vector3u");
    auto vector3f_class = bind<Vector3fD>(m, s, "Vector3f");
    auto vector3d_class = bind<Vector3dD>(m, s, "Vector3d");

    vector3f_class
        .def(py::init<const Vector3f  &>())
        .def(py::init<const Vector3fC &>())
        .def(py::init<const Vector3dD &>())
        .def(py::init<const Vector3uD &>())
        .def(py::init<const Vector3iD &>());

    vector3d_class
        .def(py::init<const Vector3d  &>())
        .def(py::init<const Vector3dC &>())
        .def(py::init<const Vector3fD &>())
        .def(py::init<const Vector3uD &>())
        .def(py::init<const Vector3iD &>());

    vector3i_class
        .def(py::init<const Vector3i  &>())
        .def(py::init<const Vector3iC &>())
        .def(py::init<const Vector3fD &>())
        .def(py::init<const Vector3dD &>())
        .def(py::init<const Vector3uD &>());

    vector3u_class
        .def(py::init<const Vector3u  &>())
        .def(py::init<const Vector3uC &>())
        .def(py::init<const Vector3fD &>())
        .def(py::init<const Vector3dD &>())
        .def(py::init<const Vector3iD &>());
}
