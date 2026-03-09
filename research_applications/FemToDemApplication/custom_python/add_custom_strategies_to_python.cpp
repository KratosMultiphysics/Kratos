//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//
// System includes


// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "custom_strategies/strategies/explicit_solver_strategy.h"
#include "custom_strategies/residualbased_DEM_coupled_newton_raphson_strategy.h"
#include "custom_strategies/femdem_residual_criteria.h"

//linear solvers
#include "linear_solvers/linear_solver.h"



namespace Kratos
{

namespace Python
{
    using namespace pybind11;

    void  AddCustomStrategiesToPython(pybind11::module& m)
    {
        typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
        typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
        typedef FemDemResidualCriteria< SparseSpaceType,  LocalSpaceType > FemDemResidualCriteriaType;
        typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
        typedef ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
        typedef BaseSolvingStrategyType::TBuilderAndSolverType BuilderAndSolverType;
        typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaType;
        typedef ResidualBasedDEMCoupledNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedDEMCoupledNewtonRaphsonStrategy;


        class_< ResidualBasedDEMCoupledNewtonRaphsonStrategy,
                typename ResidualBasedDEMCoupledNewtonRaphsonStrategy::Pointer,
                BaseSolvingStrategyType  >  (m, "ResidualBasedDEMCoupledNewtonRaphsonStrategy")
                .def(init <ModelPart& , ExplicitSolverStrategy::Pointer,  BaseSchemeType::Pointer , ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer , int  , bool , bool , bool  >())
                ;


        class_<FemDemResidualCriteria<SparseSpaceType, LocalSpaceType >,
            typename FemDemResidualCriteria<SparseSpaceType, LocalSpaceType >::Pointer,
            ConvergenceCriteriaType >
            (m,"FemDemResidualCriteria")
            .def(init< double, double>())
            ;
    }

}  // namespace Python.

} // Namespace Kratos
