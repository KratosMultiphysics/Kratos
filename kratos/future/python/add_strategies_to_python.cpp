
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
#include "containers/model.h"
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
    using BuilderType = Future::Builder<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>;
    py::class_<BuilderType, typename BuilderType::Pointer>(m, "Builder")
        .def(py::init<ModelPart &, Parameters>())
        .def("AllocateLinearSystem", py::overload_cast<LinearSystemContainer<CsrMatrix<>, SystemVector<>>&>(&BuilderType::AllocateLinearSystem))
        .def("AllocateLinearSystem", py::overload_cast<const SparseContiguousRowGraph<>&, LinearSystemContainer<CsrMatrix<>, SystemVector<>>&>(&BuilderType::AllocateLinearSystem))
        .def("AllocateLinearSystemConstraints", &BuilderType::AllocateLinearSystemConstraints)
        .def("SetUpSparseMatrixGraph", &BuilderType::SetUpSparseMatrixGraph)
        .def("SetUpMasterSlaveConstraintsGraph", &BuilderType::SetUpMasterSlaveConstraintsGraph)
        .def("ApplyLinearSystemConstraints", &BuilderType::ApplyLinearSystemConstraints)
        .def("CalculateSolutionVector", &BuilderType::CalculateSolutionVector)
        .def("Clear", &BuilderType::Clear)
        .def("GetModelPart", [&](const BuilderType &rThis) -> const ModelPart& { return rThis.GetModelPart(); }, py::return_value_policy::reference_internal)
        .def("GetEchoLevel", &BuilderType::GetEchoLevel)
        ;

    using BlockBuilderType = Future::BlockBuilder<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>;
    py::class_<BlockBuilderType, typename BlockBuilderType::Pointer>(m, "BlockBuilder")
    .def(py::init<ModelPart &, Parameters>())
    ;

    using EliminationBuilderType = Future::EliminationBuilder<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>;
    py::class_<EliminationBuilderType, typename EliminationBuilderType::Pointer>(m, "EliminationBuilder")
        .def(py::init<ModelPart &, Parameters>())
        ;

    using ImplicitSchemeType = Future::ImplicitScheme<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>;
    py::class_<ImplicitSchemeType, typename ImplicitSchemeType::Pointer>(m, "ImplicitScheme")
        .def(py::init<ModelPart &, Parameters>())
        .def("Initialize", &ImplicitSchemeType::Initialize)
        .def("InitializeSolutionStep", &ImplicitSchemeType::InitializeSolutionStep)
        .def("Predict", &ImplicitSchemeType::Predict)
        .def("FinalizeSolutionStep", &ImplicitSchemeType::FinalizeSolutionStep)
        .def("InitializeNonLinIteration", &ImplicitSchemeType::InitializeNonLinIteration)
        .def("FinalizeNonLinIteration", &ImplicitSchemeType::FinalizeNonLinIteration)
        .def("SetUpDofArrays", &ImplicitSchemeType::SetUpDofArrays)
        .def("SetUpSystemIds", &ImplicitSchemeType::SetUpSystemIds)
        // .def("Build", py::overload_cast<CsrMatrix<>, SystemVector<>>(&ImplicitSchemeType::Build)) //TODO: To be activated once we decide the parallelism
        // .def("Build", py::overload_cast<CsrMatrix<>>(&ImplicitSchemeType::Build)) //TODO: To be activated once we decide the parallelism
        // .def("Build", py::overload_cast<SystemVector<>>(&ImplicitSchemeType::Build)) //TODO: To be activated once we decide the parallelism
        .def("BuildMasterSlaveConstraints", &ImplicitSchemeType::BuildMasterSlaveConstraints)
        .def("BuildLinearSystemConstraints", &ImplicitSchemeType::BuildLinearSystemConstraints)
        .def("ApplyLinearSystemConstraints", &ImplicitSchemeType::ApplyLinearSystemConstraints)
        .def("CalculateReactions", &ImplicitSchemeType::CalculateReactions)
        .def("Update", &ImplicitSchemeType::Update)
        .def("UpdateConstraintsOnlyDofs", &ImplicitSchemeType::UpdateConstraintsOnlyDofs)
        .def("CalculateUpdateVector", &ImplicitSchemeType::CalculateUpdateVector)
        .def("Clear", &ImplicitSchemeType::Clear)
        .def("Check", &ImplicitSchemeType::Check)
        .def("GetModelPart", [&](const ImplicitSchemeType &rThis) -> const ModelPart& { return rThis.GetModelPart(); }, py::return_value_policy::reference_internal)
        .def("GetMoveMesh", &ImplicitSchemeType::GetMoveMesh)
        .def("GetReformDofsAtEachStep", &ImplicitSchemeType::GetReformDofsAtEachStep)
        .def("GetEchoLevel", &ImplicitSchemeType::GetEchoLevel)
        .def("Info", &ImplicitSchemeType::Info)
        ;

    using StaticSchemeType = Future::StaticScheme<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>;
    py::class_<StaticSchemeType, typename StaticSchemeType::Pointer, ImplicitSchemeType>(m, "StaticScheme")
        .def(py::init<ModelPart&, Parameters>())
    ;

    py::class_<Strategy, typename Strategy::Pointer>(m, "Strategy")
        .def("GetModelPart", [&](const Strategy& rThis) -> const ModelPart& {return rThis.GetModelPart();}, py::return_value_policy::reference_internal)
        .def("Info", &Strategy::Info)
        .def("Name", &Strategy::Name)
        .def("CalculateOutputData", py::overload_cast<const Variable<Vector>&>(&Strategy::CalculateOutputData, py::const_))
        .def("CalculateOutputData", py::overload_cast<const Variable<Matrix>&>(&Strategy::CalculateOutputData, py::const_))
    ;

    using LinearSolverType = Future::LinearSolver<SystemVector<>>;
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

