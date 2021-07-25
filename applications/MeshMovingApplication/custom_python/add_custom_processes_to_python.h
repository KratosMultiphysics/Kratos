//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt

#ifndef KRATOS_MESH_MOVING_PROCESSES_INCLUDED
#define KRATOS_MESH_MOVING_PROCESSES_INCLUDED

// System includes
#include "includes/define_python.h"

namespace Kratos
{
namespace Python
{

void AddCustomProcessesToPython(pybind11::module& rModule);

} // namespace Python
} // namespace Kratos

#endif // KRATOS_MESH_MOVING_PROCESSES_INCLUDED