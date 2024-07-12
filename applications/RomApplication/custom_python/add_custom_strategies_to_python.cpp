//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Raul Bravo
//
//


// System includes


// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"
#include "custom_python/add_custom_strategies_to_python.h"


#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/implicit_solving_strategy.h"
#include "custom_strategies/rom_residualbased_newton_raphson_strategy.h"
#include "custom_strategies/rom_builder_and_solver.h"
#include "custom_strategies/lspg_rom_builder_and_solver.h"
#include "custom_strategies/petrov_galerkin_rom_builder_and_solver.h"
#include "custom_strategies/global_rom_builder_and_solver.h"
#include "custom_strategies/global_petrov_galerkin_rom_builder_and_solver.h"

//linear solvers
#include "linear_solvers/linear_solver.h"


namespace Kratos {
namespace Python {

void  AddCustomStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;

    //********************************************************************
    //********************************************************************
    typedef ROMBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> ROMBuilderAndSolverType;

     py::class_<ROMBuilderAndSolverType, typename ROMBuilderAndSolverType::Pointer, BuilderAndSolverType>(m, "ROMBuilderAndSolver")
        .def(py::init< LinearSolverType::Pointer, Parameters>() )
        ;

    typedef PetrovGalerkinROMBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> PetrovGalerkinROMBuilderAndSolverType;

     py::class_<PetrovGalerkinROMBuilderAndSolverType, typename PetrovGalerkinROMBuilderAndSolverType::Pointer, ROMBuilderAndSolverType, BuilderAndSolverType>(m, "PetrovGalerkinROMBuilderAndSolver")
        .def(py::init< LinearSolverType::Pointer, Parameters>() )
        ;

    typedef GlobalROMBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> GlobalROMBuilderAndSolverType;
    typedef ResidualBasedBlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> ResidualBasedBlockBuilderAndSolverType;

    py::class_<GlobalROMBuilderAndSolverType, typename GlobalROMBuilderAndSolverType::Pointer, ResidualBasedBlockBuilderAndSolverType>(m, "GlobalROMBuilderAndSolver")
        .def(py::init< LinearSolverType::Pointer, Parameters>() )
        ;

    typedef LeastSquaresPetrovGalerkinROMBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> LeastSquaresPetrovGalerkinROMBuilderAndSolverType;

    py::class_<LeastSquaresPetrovGalerkinROMBuilderAndSolverType, typename LeastSquaresPetrovGalerkinROMBuilderAndSolverType::Pointer, GlobalROMBuilderAndSolverType>(m, "LeastSquaresPetrovGalerkinROMBuilderAndSolver")
    .def(py::init< LinearSolverType::Pointer, Parameters>() )
    .def("BuildAndApplyDirichletConditions", &LeastSquaresPetrovGalerkinROMBuilderAndSolverType::BuildAndApplyDirichletConditions)
    .def("GetRightROMBasis", &LeastSquaresPetrovGalerkinROMBuilderAndSolverType::GetRightROMBasis)
    ;

    typedef GlobalPetrovGalerkinROMBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> GlobalPetrovGalerkinROMBuilderAndSolverType;

    py::class_<GlobalPetrovGalerkinROMBuilderAndSolverType, typename GlobalPetrovGalerkinROMBuilderAndSolverType::Pointer, GlobalROMBuilderAndSolverType>(m, "GlobalPetrovGalerkinROMBuilderAndSolver")
        .def(py::init< LinearSolverType::Pointer, Parameters>() )
        ;


    //********************************************************************
    //********************************************************************

    typedef RomResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > RomResidualBasedNewtonRaphsonStrategyType;
    typedef ImplicitSolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ImplicitSolvingStrategyType;
    typedef ConvergenceCriteria<SparseSpaceType, LocalSpaceType> ConvergenceCriteriaType;

    py::class_< RomResidualBasedNewtonRaphsonStrategyType, typename RomResidualBasedNewtonRaphsonStrategyType::Pointer, ImplicitSolvingStrategyType >
        (m,"RomResidualBasedNewtonRaphsonStrategy")
        .def(py::init<ModelPart&, Parameters >() )
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, int, bool, bool, bool >())
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, ConvergenceCriteriaType::Pointer, GlobalROMBuilderAndSolverType::Pointer, int, bool, bool, bool >())
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, Parameters>())
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, ConvergenceCriteriaType::Pointer, GlobalROMBuilderAndSolverType::Pointer, Parameters>())
        .def(py::init([](ModelPart& rModelPart, BaseSchemeType::Pointer pScheme, LinearSolverType::Pointer pLinearSolver, ConvergenceCriteriaType::Pointer pConvergenceCriteria, GlobalROMBuilderAndSolverType::Pointer pBuilderAndSolver, int MaxIterations, bool CalculateReactions, bool ReformDofSetAtEachStep, bool MoveMeshFlag) {
                KRATOS_WARNING("RomResidualBasedNewtonRaphsonStrategy") << "Using deprecated constructor. Please use constructor without linear solver.";
                return std::shared_ptr<RomResidualBasedNewtonRaphsonStrategyType>(new RomResidualBasedNewtonRaphsonStrategyType(rModelPart, pScheme, pConvergenceCriteria, pBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag));
            }))
        .def(py::init([](ModelPart& rModelPart, BaseSchemeType::Pointer pScheme, LinearSolverType::Pointer pLinearSolver, ConvergenceCriteriaType::Pointer pConvergenceCriteria, GlobalROMBuilderAndSolverType::Pointer pBuilderAndSolver, Parameters Settings) {
                KRATOS_WARNING("RomResidualBasedNewtonRaphsonStrategy") << "Using deprecated constructor. Please use constructor without linear solver.";
                return std::shared_ptr<RomResidualBasedNewtonRaphsonStrategyType>(new RomResidualBasedNewtonRaphsonStrategyType(rModelPart, pScheme, pConvergenceCriteria, pBuilderAndSolver, Settings));
            }))
        .def("SetMaxIterationNumber", &RomResidualBasedNewtonRaphsonStrategyType::SetMaxIterationNumber)
        .def("GetMaxIterationNumber", &RomResidualBasedNewtonRaphsonStrategyType::GetMaxIterationNumber)
        .def("SetKeepSystemConstantDuringIterations", &RomResidualBasedNewtonRaphsonStrategyType::SetKeepSystemConstantDuringIterations)
        .def("GetKeepSystemConstantDuringIterations", &RomResidualBasedNewtonRaphsonStrategyType::GetKeepSystemConstantDuringIterations)
        .def("SetInitializePerformedFlag", &RomResidualBasedNewtonRaphsonStrategyType::SetInitializePerformedFlag)
        .def("GetInitializePerformedFlag", &RomResidualBasedNewtonRaphsonStrategyType::GetInitializePerformedFlag)
        .def("SetUseOldStiffnessInFirstIterationFlag", &RomResidualBasedNewtonRaphsonStrategyType::SetUseOldStiffnessInFirstIterationFlag)
        .def("GetUseOldStiffnessInFirstIterationFlag", &RomResidualBasedNewtonRaphsonStrategyType::GetUseOldStiffnessInFirstIterationFlag)
        .def("SetReformDofSetAtEachStepFlag", &RomResidualBasedNewtonRaphsonStrategyType::SetReformDofSetAtEachStepFlag)
        .def("GetReformDofSetAtEachStepFlag", &RomResidualBasedNewtonRaphsonStrategyType::GetReformDofSetAtEachStepFlag)
        .def("GetNonconvergedSolutions", &RomResidualBasedNewtonRaphsonStrategyType::GetNonconvergedSolutions)
        .def("SetUpNonconvergedSolutionsFlag", &RomResidualBasedNewtonRaphsonStrategyType::SetUpNonconvergedSolutionsFlag)
        ;

}

} // namespace Python.
} // Namespace Kratos

