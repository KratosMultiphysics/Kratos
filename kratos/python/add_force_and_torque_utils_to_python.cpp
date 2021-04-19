//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt

// System includes

// External includes

// Project includes
#include "python/add_force_and_torque_utils_to_python.h"
#include "utilities/force_and_torque_utils.h"

namespace Kratos
{
namespace Python
{

void AddForceAndTorqueUtilsToPython(pybind11::module& rModule)
{
    auto python_module = pybind11::class_<ForceAndTorqueUtils>(rModule, "ForceAndTorqueUtils")
        .def(pybind11::init<>())
        .def("SumForceAndTorque", &ForceAndTorqueUtils::SumForceAndTorque)
        ;
}

} /* namespace Python */
} /* namespace Kratos */