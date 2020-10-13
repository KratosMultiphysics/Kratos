//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi, Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "solving_strategies/convergence_accelerators/convergence_accelerator.h"
#include "spaces/ublas_space.h"

// Application includes
#include "custom_python/add_convergence_accelerators_to_python.h"
#include "custom_utilities/constant_relaxation_convergence_accelerator.h"
#include "custom_utilities/mvqn_convergence_accelerator.hpp"
#include "custom_utilities/mvqn_recursive_convergence_accelerator.hpp"
#include "custom_utilities/aitken_convergence_accelerator.hpp"

namespace Kratos
{

namespace Python
{

void AddConvergenceAcceleratorsToPython(pybind11::module &m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef ConvergenceAccelerator<SparseSpaceType, DenseSpaceType> BaseConvergenceAcceleratorType;

    // Constant relaxation convergence accelerator
    typedef ConstantRelaxationConvergenceAccelerator<SparseSpaceType, DenseSpaceType> ConstantRelaxationConvergenceAcceleratorType;
    typedef typename ConstantRelaxationConvergenceAcceleratorType::Pointer ConstantRelaxationConvergenceAcceleratorPointerType;
    py::class_<ConstantRelaxationConvergenceAcceleratorType, ConstantRelaxationConvergenceAcceleratorPointerType, BaseConvergenceAcceleratorType>(m, "ConstantRelaxationConvergenceAccelerator")
        .def(py::init<double>())
        .def(py::init<Parameters &>())
        .def("UpdateSolution", &ConstantRelaxationConvergenceAcceleratorType::UpdateSolution);

    // Aitken convergence accelerator
    typedef AitkenConvergenceAccelerator<SparseSpaceType, DenseSpaceType> AitkenConvergenceAcceleratorType;
    typedef typename AitkenConvergenceAcceleratorType::Pointer AitkenConvergenceAcceleratorPointerType;
    py::class_<AitkenConvergenceAcceleratorType, AitkenConvergenceAcceleratorPointerType, BaseConvergenceAcceleratorType>(m, "AitkenConvergenceAccelerator")
        .def(py::init<double>())
        .def(py::init<Parameters &>())
        .def("InitializeSolutionStep", &AitkenConvergenceAcceleratorType::InitializeSolutionStep)
        .def("UpdateSolution", &AitkenConvergenceAcceleratorType::UpdateSolution)
        .def("FinalizeNonLinearIteration", &AitkenConvergenceAcceleratorType::FinalizeNonLinearIteration)
        .def("FinalizeSolutionStep", &AitkenConvergenceAcceleratorType::FinalizeSolutionStep);

    // MVQN convergence accelerator
    typedef MVQNFullJacobianConvergenceAccelerator<SparseSpaceType, DenseSpaceType> MVQNFullJacobianConvergenceAcceleratorType;
    typedef typename MVQNFullJacobianConvergenceAcceleratorType::Pointer MVQNFullJacobianConvergenceAcceleratorPointerType;
    py::class_<MVQNFullJacobianConvergenceAcceleratorType, MVQNFullJacobianConvergenceAcceleratorPointerType, BaseConvergenceAcceleratorType>(m, "MVQNFullJacobianConvergenceAccelerator")
        .def(py::init<Parameters &>())
        .def(py::init<double, double>())
        .def("InitializeSolutionStep", &MVQNFullJacobianConvergenceAcceleratorType::InitializeSolutionStep)
        .def("UpdateSolution", &MVQNFullJacobianConvergenceAcceleratorType::UpdateSolution)
        .def("FinalizeNonLinearIteration", &MVQNFullJacobianConvergenceAcceleratorType::FinalizeNonLinearIteration)
        .def("FinalizeSolutionStep", &MVQNFullJacobianConvergenceAcceleratorType::FinalizeSolutionStep);

    // MVQN recursive convergence accelerator
    typedef MVQNRecursiveJacobianConvergenceAccelerator<SparseSpaceType, DenseSpaceType> MVQNRecursiveJacobianConvergenceAcceleratorType;
    typedef typename MVQNRecursiveJacobianConvergenceAcceleratorType::Pointer MVQNRecursiveJacobianConvergenceAcceleratorPointerType;
    py::class_<MVQNRecursiveJacobianConvergenceAcceleratorType, MVQNRecursiveJacobianConvergenceAcceleratorPointerType, BaseConvergenceAcceleratorType>(m, "MVQNRecursiveJacobianConvergenceAccelerator")
        .def(py::init<Parameters &>())
        .def(py::init<double, unsigned int, double>())
        .def("Initialize", &MVQNRecursiveJacobianConvergenceAcceleratorType::Initialize)
        .def("InitializeSolutionStep", &MVQNRecursiveJacobianConvergenceAcceleratorType::InitializeSolutionStep)
        .def("UpdateSolution", &MVQNRecursiveJacobianConvergenceAcceleratorType::UpdateSolution)
        .def("FinalizeNonLinearIteration", &MVQNRecursiveJacobianConvergenceAcceleratorType::FinalizeNonLinearIteration);
}

}  // namespace Python.

} // Namespace Kratos
