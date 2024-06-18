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
#include "custom_strategies/rom_builder_and_solver.h"
#include "custom_strategies/lspg_rom_builder_and_solver.h"
#include "custom_strategies/petrov_galerkin_rom_builder_and_solver.h"
#include "custom_strategies/global_rom_builder_and_solver.h"
#include "custom_strategies/global_petrov_galerkin_rom_builder_and_solver.h"
#include "custom_strategies/ann_prom_global_rom_builder_and_solver.h"
#include "custom_strategies/ann_prom_lspg_rom_builder_and_solver.h"
#include "custom_strategies/ann_prom_line_search_strategy.h"

/* Convergence criterias */
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

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

    typedef ConvergenceCriteria<SparseSpaceType, LocalSpaceType> ConvergenceCriteriaType;

    typedef ImplicitSolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> BaseSolvingStrategyType;

    //********************************************************************
    //********************************************************************

    typedef AnnPromLineSearchStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> AnnPromLineSearchStrategyType;

    // Line search Contact Strategy
    py::class_< AnnPromLineSearchStrategyType,
        typename AnnPromLineSearchStrategyType::Pointer,
        BaseSolvingStrategyType  >(m, "AnnPromLineSearchStrategy")
        .def(py::init < ModelPart&, Parameters >())
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, int, bool, bool, bool >())
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, Parameters >())
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, Parameters >())
        ;


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


    typedef AnnPromGlobalROMBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> AnnPromGlobalROMBuilderAndSolverType;

    py::class_<AnnPromGlobalROMBuilderAndSolverType, typename AnnPromGlobalROMBuilderAndSolverType::Pointer, ResidualBasedBlockBuilderAndSolverType>(m, "AnnPromGlobalROMBuilderAndSolver")
    .def(py::init< LinearSolverType::Pointer, Parameters>() )
    .def("SetNumberOfROMModes", &AnnPromGlobalROMBuilderAndSolverType::SetNumberOfROMModes)
    .def("SetNumberOfNNLayers", &AnnPromGlobalROMBuilderAndSolverType::SetNumberOfNNLayers)
    .def("SetNNLayer", &AnnPromGlobalROMBuilderAndSolverType::SetNNLayer)
    .def("SetPhiMatrices", &AnnPromGlobalROMBuilderAndSolverType::SetPhiMatrices)
    .def("SetReferenceSnapshot", &AnnPromGlobalROMBuilderAndSolverType::SetReferenceSnapshot)
    ;

    typedef AnnPromLeastSquaresPetrovGalerkinROMBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> AnnPromLeastSquaresPetrovGalerkinROMBuilderAndSolverType;

    py::class_<AnnPromLeastSquaresPetrovGalerkinROMBuilderAndSolverType, typename AnnPromLeastSquaresPetrovGalerkinROMBuilderAndSolverType::Pointer, ResidualBasedBlockBuilderAndSolverType>(m, "AnnPromLeastSquaresPetrovGalerkinROMBuilderAndSolver")
    .def(py::init< LinearSolverType::Pointer, Parameters>() )
    .def("SetNumberOfROMModes", &AnnPromLeastSquaresPetrovGalerkinROMBuilderAndSolverType::SetNumberOfROMModes)
    .def("SetNumberOfNNLayers", &AnnPromLeastSquaresPetrovGalerkinROMBuilderAndSolverType::SetNumberOfNNLayers)
    .def("SetNNLayer", &AnnPromLeastSquaresPetrovGalerkinROMBuilderAndSolverType::SetNNLayer)
    .def("SetPhiMatrices", &AnnPromLeastSquaresPetrovGalerkinROMBuilderAndSolverType::SetPhiMatrices)
    .def("SetReferenceSnapshot", &AnnPromLeastSquaresPetrovGalerkinROMBuilderAndSolverType::SetReferenceSnapshot)
    ;



}

} // namespace Python.
} // Namespace Kratos

