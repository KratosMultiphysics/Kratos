
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
#include "future/solving_strategies/schemes/assembly_helper.h"
#include "future/solving_strategies/schemes/implicit_scheme.h"
#include "future/solving_strategies/schemes/static_scheme.h"
#include "future/solving_strategies/strategies/implicit_strategy.h"
#include "future/solving_strategies/strategies/linear_strategy.h"

namespace Kratos::Future::Python
{

namespace py = pybind11;

void AddStrategiesToPython(py::module& m)
{
    using ImplicitAssemblyHelper = Future::AssemblyHelper<Future::ImplicitThreadLocalStorage<>, CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>;
    py::class_<ImplicitAssemblyHelper, typename ImplicitAssemblyHelper::Pointer>(m, "ImplicitAssemblyHelper")
        .def(py::init<ModelPart &, Parameters>())
        // .def("Execute",&Future::Process::Execute)
        // .def("Info",&Future::Process::Info)
        // .def("__str__", PrintObject<Future::Process>)
    ;

    using ImplicitSchemeType = Future::ImplicitScheme<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>;
    using StaticSchemeType = Future::StaticScheme<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>;
    py::class_<StaticSchemeType, typename StaticSchemeType::Pointer, ImplicitSchemeType>(m, "StaticSchemeType")
        .def(py::init<ModelPart&, Parameters>())
        // .def("Execute",&Future::Process::Execute)
        // .def("Info",&Future::Process::Info)
        // .def("__str__", PrintObject<Future::Process>)
    ;

    using LinearSolverType = Future::LinearSolver<CsrMatrix<>, SystemVector<>>;
    using LinearStrategyType = Future::LinearStrategy<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>;
    using ImplicitStrategyType = Future::ImplicitStrategy<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>;
    py::class_<LinearStrategyType, typename LinearStrategyType::Pointer, ImplicitStrategyType>(m, "LinearStrategy")
        // .def(py::init<ModelPart&, Parameters>()) //TODO: Expose this one once we fix the registry stuff
        .def(py::init<ModelPart &, typename ImplicitSchemeType::Pointer, typename LinearSolverType::Pointer, bool, bool, bool, bool>())
        .def("Initialize", &LinearStrategyType::Initialize)
        .def("InitializeSolutionStep", &LinearStrategyType::InitializeSolutionStep)
        .def("Predict", &LinearStrategyType::Predict)
        .def("SolveSolutionStep", &LinearStrategyType::SolveSolutionStep)
        .def("FinalizeSolutionStep", &LinearStrategyType::FinalizeSolutionStep)
        .def("Clear", &LinearStrategyType::Clear)
        .def("Check", &LinearStrategyType::Check)
        .def("GetEchoLevel", &LinearStrategyType::GetEchoLevel)
        .def("GetComputeReactions", &LinearStrategyType::GetComputeReactions)
        .def("GetReformDofsAtEachStep", &LinearStrategyType::GetReformDofsAtEachStep)
        .def("GetResidualNorm", &LinearStrategyType::GetResidualNorm)
        .def("GetModelPart", [&](const LinearStrategyType& rThis) -> const ModelPart& {return rThis.GetModelPart();}, py::return_value_policy::reference_internal)
        .def("SetEchoLevel", &LinearStrategyType::SetEchoLevel)
        .def("SetComputeReactions", &LinearStrategyType::SetComputeReactions)
        .def("SetReformDofsAtEachStep", &LinearStrategyType::SetReformDofsAtEachStep)
        .def("Info", &LinearStrategyType::Info);
}

}  // namespace Kratos::Future::Python.

