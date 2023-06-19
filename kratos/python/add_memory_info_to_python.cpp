//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

// External includes

// Project includes
#include "includes/memory_info.h"
#include "includes/define_python.h"
#include "add_memory_info_to_python.h"

namespace Kratos::Python
{

void AddMemoryInfoToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<MemoryInfo, MemoryInfo::Pointer>(m, "MemoryInfo")
    .def(py::init<>())
    .def_static("GetPeakMemoryUsage", &MemoryInfo::GetPeakMemoryUsage)
	.def_static("GetCurrentMemoryUsage", &MemoryInfo::GetCurrentMemoryUsage)
	.def_static("HumanReadableSize", &MemoryInfo::HumanReadableSize)
    .def("__str__", PrintObject<MemoryInfo>)
    ;
}

}  // namespace Kratos::Python.
