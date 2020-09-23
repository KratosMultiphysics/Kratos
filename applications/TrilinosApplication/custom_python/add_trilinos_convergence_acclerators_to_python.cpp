//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Philipp Bucher
//
//

// System includes

// External includes
#include "Epetra_FEVector.h"

// Project includes
#include "add_trilinos_convergence_acclerators_to_python.h"
#include "solving_strategies/convergence_accelerators/convergence_accelerator.h"
#include "spaces/ublas_space.h"
#include "trilinos_space.h"

namespace Kratos {
namespace Python {

void AddTrilinosConvergenceAcceleratorsToPython(pybind11::module &m)
{
    namespace py = pybind11;

    using TrilinosSparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
    using TrilinosLocalSpaceType  = UblasSpace<double, Matrix, Vector>;

    // Convergence accelerator base class
    using TrilinosBaseConvergenceAcceleratorType = ConvergenceAccelerator<TrilinosSparseSpaceType, TrilinosLocalSpaceType>;

    py::class_<TrilinosBaseConvergenceAcceleratorType, typename TrilinosBaseConvergenceAcceleratorType::Pointer>(m, "TrilinosConvergenceAccelerator")
        .def(py::init<Parameters>())
        .def("Initialize",                   &TrilinosBaseConvergenceAcceleratorType::Initialize)
        .def("InitializeSolutionStep",       &TrilinosBaseConvergenceAcceleratorType::InitializeSolutionStep)
        .def("InitializeNonLinearIteration", &TrilinosBaseConvergenceAcceleratorType::InitializeNonLinearIteration)
        .def("UpdateSolution",               &TrilinosBaseConvergenceAcceleratorType::UpdateSolution)
        .def("FinalizeNonLinearIteration",   &TrilinosBaseConvergenceAcceleratorType::FinalizeNonLinearIteration)
        .def("FinalizeSolutionStep",         &TrilinosBaseConvergenceAcceleratorType::FinalizeSolutionStep)
        .def("Finalize",                     &TrilinosBaseConvergenceAcceleratorType::Finalize)
        .def("GetEchoLevel",                 &TrilinosBaseConvergenceAcceleratorType::GetEchoLevel)
        .def("SetEchoLevel",                 &TrilinosBaseConvergenceAcceleratorType::SetEchoLevel)
    ;

}

}  // namespace Python.
} // Namespace Kratos
