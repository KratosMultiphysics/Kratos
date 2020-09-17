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
//
//

// System includes

// External includes

// Project includes
#include "add_convergence_accelerators_to_python.h"
#include "includes/define.h"
#include "solving_strategies/convergence_accelerators/convergence_accelerator.h"
#include "spaces/ublas_space.h"

namespace Kratos
{

namespace Python
{

void AddConvergenceAcceleratorsToPython(pybind11::module &m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, Matrix, Vector > TSpace;
    typedef ConvergenceAccelerator< TSpace > BaseConvergenceAcceleratorType;

    // Convergence accelerator base class
    py::class_<ConvergenceAccelerator<TSpace>>(m, "ConvergenceAccelerator")
        .def(py::init<>())
        .def("Initialize", &ConvergenceAccelerator<TSpace>::Initialize)
        .def("InitializeSolutionStep", &ConvergenceAccelerator<TSpace>::InitializeSolutionStep)
        .def("InitializeNonLinearIteration", &ConvergenceAccelerator<TSpace>::InitializeNonLinearIteration)
        .def("UpdateSolution", &ConvergenceAccelerator<TSpace>::UpdateSolution)
        .def("FinalizeNonLinearIteration", &ConvergenceAccelerator<TSpace>::FinalizeNonLinearIteration)
        .def("FinalizeSolutionStep", &ConvergenceAccelerator<TSpace>::FinalizeSolutionStep)
        .def("SetEchoLevel", &ConvergenceAccelerator<TSpace>::SetEchoLevel);

}

}  // namespace Python.

} // Namespace Kratos
