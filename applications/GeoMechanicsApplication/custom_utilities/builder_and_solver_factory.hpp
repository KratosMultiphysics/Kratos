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
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace, class TLinearSolver>
class BuilderAndSolverFactory
{

public:
    using BuilderAndSolverType = BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>;

    static std::shared_ptr<BuilderAndSolverType> Create(const Parameters& rSolverSettings,
                                                        typename TLinearSolver::Pointer pNewLinearSystemSolver)
    {
        const std::string block_builder_entry = "block_builder";
        KRATOS_ERROR_IF_NOT(rSolverSettings.Has(block_builder_entry)) <<
        "the block_builder parameter is not defined, aborting BuilderAndSolverCreation";

        if (rSolverSettings[block_builder_entry].GetBool())
        {
            return std::make_shared<ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>>(pNewLinearSystemSolver);
        }

        return std::make_shared<ResidualBasedEliminationBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>>(pNewLinearSystemSolver);
    }
};


}
