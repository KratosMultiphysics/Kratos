//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt

#ifndef KRATOS_ADD_FORCE_AND_TORQUE_UTILS_TO_PYTHON
#define KRATOS_ADD_FORCE_AND_TORQUE_UTILS_TO_PYTHON

// System includes
#include <pybind11/pybind11.h>

// External includes


// Project includes

namespace Kratos
{
namespace Python
{

void AddForceAndTorqueUtilsToPython(pybind11::module& rModule);

} /* namespace Python */
} /* namespace Kratos */

#endif // KRATOS_ADD_FORCE_AND_TORQUE_UTILS_TO_PYTHON defined