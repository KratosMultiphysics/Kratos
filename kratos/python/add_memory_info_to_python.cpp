//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/memory_info.h"

namespace Kratos
{

namespace Python
{

//
void  AddMemoryInfoToPython(pybind11::module& m)
{
    using namespace pybind11;

    class_<MemoryInfo, MemoryInfo::Pointer>(m, "MemoryInfo")
    .def(init<>())
    .def_static("GetPeakMemoryUsage", &MemoryInfo::GetPeakMemoryUsage)
	.def_static("GetCurrentMemoryUsage", &MemoryInfo::GetCurrentMemoryUsage)
    ;
}

}  // namespace Python.

} // Namespace Kratos
