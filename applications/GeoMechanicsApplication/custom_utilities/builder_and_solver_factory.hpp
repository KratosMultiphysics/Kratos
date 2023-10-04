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
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace, class TLinearSolver>
class BuilderAndSolverFactory
{

public:
    static std::shared_ptr<BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>> Create(const Parameters& rSolverSettings, typename TLinearSolver::Pointer pNewLinearSystemSolver)
    {
        if (!rSolverSettings.Has("block_builder"))
        {
            return nullptr;
        }

        if (rSolverSettings["block_builder"].GetBool())
        {
            return std::make_shared<ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>>(pNewLinearSystemSolver);
        }


        return nullptr;
    }
};


}
