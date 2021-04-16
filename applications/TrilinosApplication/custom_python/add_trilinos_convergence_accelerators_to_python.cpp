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
#include "solving_strategies/convergence_accelerators/convergence_accelerator.h"
#include "spaces/ublas_space.h"
#include "trilinos_space.h"

// Application includes
#include "trilinos_pointer_wrapper.h"
#include "add_trilinos_convergence_accelerators_to_python.h"
#include "custom_utilities/trilinos_mvqn_recursive_convergence_accelerator.hpp"
#include "../FSIApplication/custom_utilities/aitken_convergence_accelerator.hpp" // To be removed

namespace Kratos {
namespace Python {

typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;
typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;

void AuxiliarUpdateSolution(
    ConvergenceAccelerator<TrilinosSparseSpaceType, TrilinosLocalSpaceType> &dummy,
    AuxiliaryVectorWrapper &rResidualVector,
    AuxiliaryVectorWrapper &rIterationGuess)
{
    dummy.UpdateSolution(rResidualVector.GetReference(), rIterationGuess.GetReference());
}

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
        .def("UpdateSolution",               AuxiliarUpdateSolution)
        .def("FinalizeNonLinearIteration",   &TrilinosBaseConvergenceAcceleratorType::FinalizeNonLinearIteration)
        .def("FinalizeSolutionStep",         &TrilinosBaseConvergenceAcceleratorType::FinalizeSolutionStep)
        .def("Finalize",                     &TrilinosBaseConvergenceAcceleratorType::Finalize)
        .def("GetEchoLevel",                 &TrilinosBaseConvergenceAcceleratorType::GetEchoLevel)
        .def("SetEchoLevel",                 &TrilinosBaseConvergenceAcceleratorType::SetEchoLevel)
    ;

    // Convergence accelerators (from FSIApplication)
    typedef AitkenConvergenceAccelerator<TrilinosSparseSpaceType, TrilinosLocalSpaceType> TrilinosAitkenAccelerator;
    typedef typename TrilinosAitkenAccelerator::Pointer TrilinosAitkenAcceleratorPointer;
    py::class_<TrilinosAitkenAccelerator, TrilinosAitkenAcceleratorPointer, TrilinosBaseConvergenceAcceleratorType>(m,"TrilinosAitkenConvergenceAccelerator")
        .def(py::init<double>())
        .def(py::init< Parameters>())
        ;

    typedef TrilinosMVQNRecursiveJacobianConvergenceAccelerator<TrilinosSparseSpaceType, TrilinosLocalSpaceType> TrilinosMVQNRecursiveAccelerator;
    typedef typename TrilinosMVQNRecursiveAccelerator::Pointer TrilinosMVQNRecursiveAcceleratorPointer;
    py::class_< TrilinosMVQNRecursiveAccelerator, TrilinosMVQNRecursiveAcceleratorPointer, TrilinosBaseConvergenceAcceleratorType>(m,"TrilinosMVQNRecursiveJacobianConvergenceAccelerator")
        .def(py::init< ModelPart&, const Epetra_MpiComm&, Parameters >())
        .def(py::init< ModelPart&, const Epetra_MpiComm&, double, unsigned int >())
        ;

}

} // namespace Python.
} // Namespace Kratos
