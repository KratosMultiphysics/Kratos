//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Klaus B. Sautter
//


// Project includes
#include "custom_python/add_custom_processes_to_python.h"

//Processes
#include "custom_processes/sliding_edge_process.h"
#include "custom_processes/edge_cable_element_process.h"
#include "custom_processes/apply_weak_sliding_process.h"
#include "custom_processes/empirical_spring_element_process.h"

namespace Kratos {
namespace Python {

void  AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    /// Processes
    py::class_<SlidingEdgeProcess, SlidingEdgeProcess::Pointer, Process>(m,"SlidingEdgeProcess")
        .def(py::init<ModelPart&,Parameters>())
        ;

    py::class_<ApplyWeakSlidingProcess, ApplyWeakSlidingProcess::Pointer, Process>(m,"ApplyWeakSlidingProcess")
        .def(py::init<ModelPart&,Parameters>())
        ;

    py::class_<EdgeCableElementProcess, EdgeCableElementProcess::Pointer, Process>(m,"EdgeCableElementProcess")
        .def(py::init<ModelPart&,Parameters>())
        ;

    py::class_<EmpiricalSpringElementProcess, EmpiricalSpringElementProcess::Pointer, Process>(m,"EmpiricalSpringElementProcess")
        .def(py::init<ModelPart&,Parameters,std::vector<double>>())
        ;
}

}  // namespace Python.
} // Namespace Kratos

