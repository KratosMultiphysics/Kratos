// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Manuel Messmer
//

// System includes

// External includes


// Project includes
#include "custom_python/add_custom_processes_to_python.h"

//Processes
#include "custom_processes/perturb_geometry_subgrid_process.h"
#include "custom_processes/perturb_geometry_sparse_process.h"

namespace Kratos {
namespace Python {

void AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    /// Processes
    py::class_<PerturbGeometrySubgridProcess, PerturbGeometrySubgridProcess::Pointer, Process>(m,"PerturbGeometrySubgridProcess")
        .def(py::init<ModelPart&,double>())
        .def("SetEchoLevel", &PerturbGeometrySubgridProcess::SetEchoLevel)
        .def("CreateEigenvectors", &PerturbGeometrySubgridProcess::CreateEigenvectors)
        .def("AssembleEigenvectors", &PerturbGeometrySubgridProcess::AssembleEigenvectors)
        ;

    py::class_<PerturbGeometrySparseProcess, PerturbGeometrySparseProcess::Pointer, Process>(m,"PerturbGeometrySparseProcess")
        .def(py::init<ModelPart&,double>())
        .def("SetEchoLevel", &PerturbGeometrySparseProcess::SetEchoLevel)
        .def("CreateEigenvectors", &PerturbGeometrySparseProcess::CreateEigenvectors)
        .def("AssembleEigenvectors", &PerturbGeometrySparseProcess::AssembleEigenvectors)
        ;
}

}  // namespace Python.
} // Namespace Kratos

