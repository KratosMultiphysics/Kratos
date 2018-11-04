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
#include "spaces/ublas_space.h"

#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_convergence_accelerators_to_python.h"
#include "custom_utilities/convergence_accelerator.hpp"
#include "custom_utilities/mvqn_convergence_accelerator.hpp"
#include "custom_utilities/mvqn_recursive_convergence_accelerator.hpp"
#include "custom_utilities/aitken_convergence_accelerator.hpp"

namespace Kratos
{

namespace Python
{

void AddConvergenceAcceleratorsToPython(pybind11::module &m)
{
    using namespace pybind11;

    typedef UblasSpace<double, Matrix, Vector > TSpace;
    typedef ConvergenceAccelerator< TSpace > BaseConvergenceAcceleratorType;
    
    // Convergence accelerator base class
    class_<ConvergenceAccelerator<TSpace>>(m, "ConvergenceAccelerator")
        .def(init<>())
        .def("Initialize", &ConvergenceAccelerator<TSpace>::Initialize)
        .def("InitializeSolutionStep", &ConvergenceAccelerator<TSpace>::InitializeSolutionStep)
        .def("InitializeNonLinearIteration", &ConvergenceAccelerator<TSpace>::InitializeNonLinearIteration)
        .def("UpdateSolution", &ConvergenceAccelerator<TSpace>::UpdateSolution)
        .def("FinalizeNonLinearIteration", &ConvergenceAccelerator<TSpace>::FinalizeNonLinearIteration)
        .def("FinalizeSolutionStep", &ConvergenceAccelerator<TSpace>::FinalizeSolutionStep)
        .def("SetEchoLevel", &ConvergenceAccelerator<TSpace>::SetEchoLevel);

    // Aitken convergence accelerator
    class_<AitkenConvergenceAccelerator<TSpace>, BaseConvergenceAcceleratorType>(m, "AitkenConvergenceAccelerator")
        .def(init<double>())
        .def(init<Parameters &>())
        .def("InitializeSolutionStep", &AitkenConvergenceAccelerator<TSpace>::InitializeSolutionStep)
        .def("UpdateSolution", &AitkenConvergenceAccelerator<TSpace>::UpdateSolution)
        .def("FinalizeNonLinearIteration", &AitkenConvergenceAccelerator<TSpace>::FinalizeNonLinearIteration)
        .def("FinalizeSolutionStep", &AitkenConvergenceAccelerator<TSpace>::FinalizeSolutionStep);

    // MVQN convergence accelerator
    class_<MVQNFullJacobianConvergenceAccelerator<TSpace>, BaseConvergenceAcceleratorType>(m, "MVQNFullJacobianConvergenceAccelerator")
        .def(init<double>())
        .def(init<Parameters &>())
        .def("InitializeSolutionStep", &MVQNFullJacobianConvergenceAccelerator<TSpace>::InitializeSolutionStep)
        .def("UpdateSolution", &MVQNFullJacobianConvergenceAccelerator<TSpace>::UpdateSolution)
        .def("FinalizeNonLinearIteration", &MVQNFullJacobianConvergenceAccelerator<TSpace>::FinalizeNonLinearIteration)
        .def("FinalizeSolutionStep", &MVQNFullJacobianConvergenceAccelerator<TSpace>::FinalizeSolutionStep);

    // MVQN recursive convergence accelerator
    class_<MVQNRecursiveJacobianConvergenceAccelerator<TSpace>, BaseConvergenceAcceleratorType>(m, "MVQNRecursiveJacobianConvergenceAccelerator")
        .def(init<double, unsigned int>())
        .def(init<Parameters &>())
        .def(init<double, unsigned int, double>())
        .def("Initialize", &MVQNRecursiveJacobianConvergenceAccelerator<TSpace>::Initialize)
        .def("InitializeSolutionStep", &MVQNRecursiveJacobianConvergenceAccelerator<TSpace>::InitializeSolutionStep)
        .def("UpdateSolution", &MVQNRecursiveJacobianConvergenceAccelerator<TSpace>::UpdateSolution)
        .def("FinalizeNonLinearIteration", &MVQNRecursiveJacobianConvergenceAccelerator<TSpace>::FinalizeNonLinearIteration);
}

}  // namespace Python.

} // Namespace Kratos

