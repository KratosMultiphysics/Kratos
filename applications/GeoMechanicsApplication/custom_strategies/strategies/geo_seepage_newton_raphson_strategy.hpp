// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//

#pragma once

#include <string>
#include <vector>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

// A Newton-Raphson strategy that can switch seepage boundary conditions between Dirichlet and
// Neumann while iterating.
//
// Applying a switch is the condition's job: the scheme calls
// Condition::InitializeNonLinearIteration on every non-linear iteration, and GeoSeepageCondition
// fixes or frees its nodes there. This strategy owns the two things the condition cannot do:
// deciding when to switch, and making sure the solver does not declare convergence on an iteration
// whose boundary configuration has just changed underneath it.
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class GeoSeepageNewtonRaphsonStrategy
    : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(GeoSeepageNewtonRaphsonStrategy);

    using BaseType   = ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;
    using MotherType = ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;
    using TConvergenceCriteriaType = ConvergenceCriteria<TSparseSpace, TDenseSpace>;
    using TBuilderAndSolverType    = typename BaseType::TBuilderAndSolverType;
    using TSchemeType              = typename BaseType::TSchemeType;
    using DofsArrayType            = typename BaseType::DofsArrayType;
    using TSystemMatrixType        = typename BaseType::TSystemMatrixType;
    using TSystemVectorType        = typename BaseType::TSystemVectorType;

    // The base class is a dependent template base, so every inherited data member used by the
    // copied SolveSolutionStep must be pulled into scope explicitly.
    using MotherType::mCalculateReactionsFlag;
    using MotherType::mMaxIterationNumber;
    using MotherType::mNonconvergedSolutionsMatrix;
    using MotherType::mpA;  // Tangent matrix
    using MotherType::mpb;  // Residual vector of iteration i
    using MotherType::mpConvergenceCriteria;
    using MotherType::mpDx; // Delta x of iteration i
    using MotherType::mStoreNonconvergedSolutionsFlag;
    using MotherType::mUseOldStiffnessInFirstIteration;

    GeoSeepageNewtonRaphsonStrategy(ModelPart&                                 rModelPart,
                                    typename TSchemeType::Pointer              pScheme,
                                    typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
                                    typename TBuilderAndSolverType::Pointer    pNewBuilderAndSolver,
                                    int  MaxIterations          = 30,
                                    bool CalculateReactions     = false,
                                    bool ReformDofSetAtEachStep = false,
                                    bool MoveMeshFlag           = false)
        : MotherType(rModelPart,
                     pScheme,
                     pNewConvergenceCriteria,
                     pNewBuilderAndSolver,
                     MaxIterations,
                     CalculateReactions,
                     ReformDofSetAtEachStep,
                     MoveMeshFlag)
    {
    }

    [[nodiscard]] std::string Info() const override { return "GeoSeepageNewtonRaphsonStrategy"; }

protected:
    // Decides whether any seepage condition must change its boundary type, and applies that
    // decision to the conditions' Properties. Returns true if at least one condition changed.
    //
    // Step 3 of the prototype implements this. Until then nothing ever switches, which makes this
    // strategy behave exactly like its base class.
    virtual bool UpdateSeepageBoundaryConditions() { return false; }

    // Re-establishes the degree of freedom set and the system vectors after a switch changed which
    // degrees of freedom are free.
    //
    // Step 4 of the prototype implements this. Note that under a block builder and solver this may
    // turn out to be unnecessary, because that builder gives every degree of freedom an equation id
    // and re-applies the Dirichlet conditions on every build.
    virtual void RebuildSystem() {}
}; // Class GeoSeepageNewtonRaphsonStrategy

} // namespace Kratos

