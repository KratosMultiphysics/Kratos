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

// Project includes
#include "add_convergence_accelerators_to_python.h"
#include "solving_strategies/convergence_accelerators/convergence_accelerator.h"
#include "spaces/ublas_space.h"

namespace Kratos::Python {

void AddConvergenceAcceleratorsToPython(pybind11::module &m)
{
    namespace py = pybind11;

    using SparseSpaceType = UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>>;
    using LocalSpaceType  = UblasSpace<double, Matrix, Vector>;

    // Convergence accelerator base class
    using BaseConvergenceAcceleratorType = ConvergenceAccelerator<SparseSpaceType, LocalSpaceType>;

    py::class_<BaseConvergenceAcceleratorType, typename BaseConvergenceAcceleratorType::Pointer>(m, "ConvergenceAccelerator")
        .def(py::init<Parameters>())
        .def("Initialize",                   &BaseConvergenceAcceleratorType::Initialize)
        .def("InitializeSolutionStep",       &BaseConvergenceAcceleratorType::InitializeSolutionStep)
        .def("InitializeNonLinearIteration", &BaseConvergenceAcceleratorType::InitializeNonLinearIteration)
        .def("UpdateSolution",               &BaseConvergenceAcceleratorType::UpdateSolution)
        .def("FinalizeNonLinearIteration",   &BaseConvergenceAcceleratorType::FinalizeNonLinearIteration)
        .def("FinalizeSolutionStep",         &BaseConvergenceAcceleratorType::FinalizeSolutionStep)
        .def("Finalize",                     &BaseConvergenceAcceleratorType::Finalize)
        .def("GetEchoLevel",                 &BaseConvergenceAcceleratorType::GetEchoLevel)
        .def("SetEchoLevel",                 &BaseConvergenceAcceleratorType::SetEchoLevel)
        .def("IsBlockNewton",                &BaseConvergenceAcceleratorType::IsBlockNewton) 
    ;

}

}  // namespace Kratos::Python.
