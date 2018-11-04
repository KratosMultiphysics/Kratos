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



// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/process_info.h"
#include "python/add_process_info_to_python.h"
#include "containers/data_value_container.h"

namespace Kratos
{

namespace Python
{


ProcessInfo::Pointer ProcessInfoGetPreviousSolutionStepInfo(ProcessInfo & rProcessInfo)
{
	return rProcessInfo.pGetPreviousSolutionStepInfo();
}

//
void  AddProcessInfoToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<ProcessInfo, ProcessInfo::Pointer, DataValueContainer, Flags >(m,"ProcessInfo")
    .def(py::init<>())
    .def("CreateSolutionStepInfo", &ProcessInfo::CreateSolutionStepInfo)
	.def("GetPreviousSolutionStepInfo", ProcessInfoGetPreviousSolutionStepInfo)
    .def("__str__", PrintObject<ProcessInfo>)
    ;
}
}  // namespace Python.

} // Namespace Kratos

