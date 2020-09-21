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

#if !defined(KRATOS_CONVERGENCE_ACCELERATOR_PYTHON_H_INCLUDED )
#define  KRATOS_CONVERGENCE_ACCELERATOR_PYTHON_H_INCLUDED

// System includes
#include <pybind11/pybind11.h>

// External includes


// Project includes
#include "solving_strategies/convergence_accelerators/convergence_accelerator.h"


namespace Kratos {
namespace Python {

// function to expose the base convergence accelerator to python, depending on the space type
template<class TSparseSpace, class TDenseSpace>
void AddBaseConvergenceAcceleratorToPython(pybind11::module& m, const std::string& rName)
{
    namespace py = pybind11;

    using BaseConvergenceAcceleratorType = ConvergenceAccelerator<TSparseSpace, TDenseSpace>;

    py::class_<BaseConvergenceAcceleratorType, typename BaseConvergenceAcceleratorType::Pointer>(m, rName.c_str())
        .def(py::init<>())
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
    ;
}

// function to expose a convergence accelerator to python, depending on the space type
template<class TSparseSpace, class TDenseSpace, template<class SparseSpace, class DenseSpace> class TConvergenceAccelerator >
void AddConvergenceAcceleratorToPython(pybind11::module& m, const std::string& rName)
{
    namespace py = pybind11;

    using BaseConvergenceAcceleratorType = ConvergenceAccelerator<TSparseSpace, TDenseSpace>;
    using ConvergenceAcceleratorType     = TConvergenceAccelerator<TSparseSpace, TDenseSpace>;

    py::class_<ConvergenceAcceleratorType, typename ConvergenceAcceleratorType::Pointer, BaseConvergenceAcceleratorType>(m, rName.c_str())
        .def(py::init<>())
        .def(py::init<Parameters>())
    ;
}

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_CONVERGENCE_ACCELERATOR_PYTHON_H_INCLUDED  defined
