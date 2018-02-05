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
#include <boost/python.hpp>

// Project includes
#include "includes/memory_info.h"
using namespace boost::python;

namespace Kratos
{

namespace Python
{

//
void  AddMemoryInfoToPython()
{
    using namespace boost::python;

    class_<MemoryInfo, MemoryInfo::Pointer, boost::noncopyable>("MemoryInfo")
    .def(init<>())
    .def("GetPeakMemoryUsage", &MemoryInfo::GetPeakMemoryUsage)
    .staticmethod("GetPeakMemoryUsage")
	.def("GetCurrentMemoryUsage", &MemoryInfo::GetCurrentMemoryUsage)
    .staticmethod("GetCurrentMemoryUsage")
    .def(self_ns::str(self))
    ;
}

}  // namespace Python.

} // Namespace Kratos
