/*
//  KRATOS _______
//        / ____(_)___ ____  ____
//       / __/ / / __ `/ _ \/ __ \
//      / /___/ / /_/ /  __/ / / /
//     /_____/_/\__, /\___/_/ /_/ SolversApplication
//             /____/
//
//  Author: Manuel Messmer
*/

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

