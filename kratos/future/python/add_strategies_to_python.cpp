
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
#include "future/solving_strategies/builders/builder.h"
#include "future/solving_strategies/schemes/implicit_scheme.h"
#include "future/solving_strategies/schemes/static_scheme.h"
#include "future/solving_strategies/strategies/strategy.h"
#include "future/solving_strategies/strategies/implicit_strategy.h"
#include "future/solving_strategies/strategies/linear_strategy.h"

namespace Kratos::Future::Python
{

namespace py = pybind11;

void AddStrategiesToPython(py::module& m)
{
    using ImplicitBuilder = Future::Builder<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>;
    py::class_<ImplicitBuilder, typename ImplicitBuilder::Pointer>(m, "ImplicitBuilder")
        .def(py::init<ModelPart &, Parameters>())
        // .def("Execute",&Future::Process::Execute)
        // .def("Info",&Future::Process::Info)
        // .def("__str__", PrintObject<Future::Process>)
    ;

    using ImplicitSchemeType = Future::ImplicitScheme<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>;
    py::class_<ImplicitSchemeType, typename ImplicitSchemeType::Pointer>(m, "ImplicitScheme")
        .def(py::init<ModelPart&, Parameters>())
        // .def("Execute",&Future::Process::Execute)
        // .def("Info",&Future::Process::Info)
        // .def("__str__", PrintObject<Future::Process>)
    ;

    using StaticSchemeType = Future::StaticScheme<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>;
    py::class_<StaticSchemeType, typename StaticSchemeType::Pointer, ImplicitSchemeType>(m, "StaticScheme")
        .def(py::init<ModelPart&, Parameters>())
        // .def("Execute",&Future::Process::Execute)
        // .def("Info",&Future::Process::Info)
        // .def("__str__", PrintObject<Future::Process>)
    ;

    py::class_<Strategy, typename Strategy::Pointer>(m, "Strategy")
        .def("GetModelPart", [&](const Strategy& rThis) -> const ModelPart& {return rThis.GetModelPart();}, py::return_value_policy::reference_internal)
        .def("Info", &Strategy::Info)
        .def("Name", &Strategy::Name)
    ;

    using LinearSolverType = Future::LinearSolver<CsrMatrix<>, SystemVector<>>;
    using ImplicitStrategyType = Future::ImplicitStrategy<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>;
    py::class_<ImplicitStrategyType, typename ImplicitStrategyType::Pointer, Strategy>(m, "ImplicitStrategy")
        // .def(py::init<ModelPart&, Parameters>()) //TODO: Expose this one once we fix the registry stuff
        .def(py::init<ModelPart &, typename ImplicitSchemeType::Pointer, typename LinearSolverType::Pointer, bool, bool, bool, bool>())
        .def("Initialize", &ImplicitStrategyType::Initialize)
        .def("InitializeSolutionStep", &ImplicitStrategyType::InitializeSolutionStep)
        .def("Predict", &ImplicitStrategyType::Predict)
        .def("SolveSolutionStep", &ImplicitStrategyType::SolveSolutionStep)
        .def("FinalizeSolutionStep", &ImplicitStrategyType::FinalizeSolutionStep)
        .def("Clear", &ImplicitStrategyType::Clear)
        .def("Check", &ImplicitStrategyType::Check)
        .def("GetEchoLevel", &ImplicitStrategyType::GetEchoLevel)
        .def("GetComputeReactions", &ImplicitStrategyType::GetComputeReactions)
        .def("GetReformDofsAtEachStep", &ImplicitStrategyType::GetReformDofsAtEachStep)
        .def("GetResidualNorm", &ImplicitStrategyType::GetResidualNorm)
        .def("SetEchoLevel", &ImplicitStrategyType::SetEchoLevel)
        .def("SetComputeReactions", &ImplicitStrategyType::SetComputeReactions)
        .def("SetReformDofsAtEachStep", &ImplicitStrategyType::SetReformDofsAtEachStep)
    ;

    using LinearStrategyType = Future::LinearStrategy<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>;
    py::class_<LinearStrategyType, typename LinearStrategyType::Pointer, ImplicitStrategyType>(m, "LinearStrategy")
        // .def(py::init<ModelPart&, Parameters>()) //TODO: Expose this one once we fix the registry stuff
        .def(py::init<ModelPart &, typename ImplicitSchemeType::Pointer, typename LinearSolverType::Pointer, bool, bool, bool, bool>())
    ;
}

}  // namespace Kratos::Future::Python.

