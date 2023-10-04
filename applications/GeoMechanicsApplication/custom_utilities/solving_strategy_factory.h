// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#pragma once

#include <memory>
#include <string>

#include "solving_strategies/strategies/solving_strategy.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/convergencecriterias/mixed_generic_criteria.h"

namespace Kratos
{

class SolvingStrategyFactory
{
public:
    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
    using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;
    using MixedGenericCriteriaType = MixedGenericCriteria<SparseSpaceType, LocalSpaceType>;
    using ConvergenceVariableListType =  MixedGenericCriteriaType::ConvergenceVariableListType;
    using ConvergenceCriteriaType = ConvergenceCriteria<SparseSpaceType, LocalSpaceType>;

    static [[nodiscard]] std::unique_ptr<SolvingStrategy<SparseSpaceType, LocalSpaceType>> Create(Parameters& rSolverSettings, ModelPart& rModelPart);
};

}
