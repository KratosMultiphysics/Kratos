//
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Angel Celigueta
//


// System includes

// External includes

// Project includes
#include "includes/define_python.h"

#include "add_quaternion_to_python.h"
#include "utilities/quaternion.h"

namespace Kratos {
namespace Python {
/**
 * @brief Sets the current table utility on the process info
 * @param rCurrentProcessInfo The process info
 */
namespace py = pybind11;

double QuaternionGetX(Quaternion<double>& ThisQuaternion) { return ThisQuaternion.X(); }
double QuaternionGetY(Quaternion<double>& ThisQuaternion) { return ThisQuaternion.Y(); }
double QuaternionGetZ(Quaternion<double>& ThisQuaternion) { return ThisQuaternion.Z(); }
double QuaternionGetW(Quaternion<double>& ThisQuaternion) { return ThisQuaternion.W(); }
void QuaternionSetX(Quaternion<double>& ThisQuaternion, double Value) { ThisQuaternion.SetX(Value); }
void QuaternionSetY(Quaternion<double>& ThisQuaternion, double Value) { ThisQuaternion.SetY(Value); }
void QuaternionSetZ(Quaternion<double>& ThisQuaternion, double Value) { ThisQuaternion.SetZ(Value); }
void QuaternionSetW(Quaternion<double>& ThisQuaternion, double Value) { ThisQuaternion.SetW(Value); }

void  AddQuaternionToPython(pybind11::module& m) {

    py::class_<Quaternion<double>, Quaternion<double>::Pointer>(m, "Quaternion")
    .def(py::init<>())
    .def("__str__", PrintObject<Quaternion<double>>)
    .def_property("X", QuaternionGetX, QuaternionSetX)
    .def_property("Y", QuaternionGetY, QuaternionSetY)
    .def_property("Z", QuaternionGetZ, QuaternionSetZ)
    .def_property("W", QuaternionGetW, QuaternionSetW)
    .def("Normalize",        &Quaternion<double>::normalize)
    .def_static("Identity",  &Quaternion<double>::Identity)
    .def("ToRotationMatrix", &Quaternion<double>::ToRotationMatrix<Matrix>)
    .def("MultiplyQuaternion", [](Quaternion<double>& self, Quaternion<double>& a, Quaternion<double>& b){return Quaternion<double>(operator*(a,b));})
    ;
}

} // namespace Python.
} // Namespace Kratos
