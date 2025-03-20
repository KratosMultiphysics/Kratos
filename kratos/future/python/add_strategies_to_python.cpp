
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
#include "includes/define_python.h"

// Future Extensions
#include "future/python/add_strategies_to_python.h"
#include "future/solving_strategies/schemes/new_scheme.h"
#include "future/solving_strategies/strategies/linear_strategy.h"

namespace Kratos::Future::Python
{

namespace py = pybind11;

void AddStrategiesToPython(py::module& m)
{
    using SchemeType = Future::NewScheme<CsrMatrix<>, SystemVector<>>;
    py::class_<SchemeType, typename SchemeType::Pointer, Flags>(m,"NewScheme")
        // .def(py::init<>())
        // .def(py::init<ModelPart&, Parameters>())
        // .def("Execute",&Future::Process::Execute)
        // .def("Info",&Future::Process::Info)
        // .def("__str__", PrintObject<Future::Process>)
    ;

    // Pack everything under the "Future" submodule
    using LinearSolverType = Future::LinearSolver<CsrMatrix<>, SystemVector<>>;
    using LinearStrategyType = Future::LinearStrategy<CsrMatrix<>, SystemVector<>, LinearSolverType>;
    py::class_<LinearStrategyType, typename LinearStrategyType::Pointer, Flags>(m,"LinearStrategy")
        // .def(py::init<>())
        // .def(py::init<ModelPart&, Parameters>())
        // .def("Execute",&Future::Process::Execute)
        // .def("Info",&Future::Process::Info)
        // .def("__str__", PrintObject<Future::Process>)
    ;
}

}  // namespace Kratos::Future::Python.

