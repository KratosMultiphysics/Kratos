//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez
//

#if !defined(KRATOS_UPWIND_RESIDUALBASED_NEWTON_RAPHSON_STRATEGY)
#define  KRATOS_UPWIND_RESIDUALBASED_NEWTON_RAPHSON_STRATEGY
// System includes

// External includes

// Project includes
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

namespace Kratos
{
template<class TSparseSpace,
class TDenseSpace, // = DenseSpace<double>,
class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
>
class UpwindResidualBasedNewtonRaphsonStrategy
: public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(UpwindResidualBasedNewtonRaphsonStrategy);

    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef typename BaseType::TSchemeType TSchemeType;

    ///@}
    ///@name Life Cycle

    ///@{

    /**
    * Constructor.
    */

    // constructor with Builder and Solver
    UpwindResidualBasedNewtonRaphsonStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
    )
    : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme,
        pNewLinearSolver,pNewConvergenceCriteria,MaxIterations,CalculateReactions,ReformDofSetAtEachStep,
        MoveMeshFlag)
    {
    }

    // Destructor
    ~UpwindResidualBasedNewtonRaphsonStrategy() = default;

private:


}; /* Class UpwindResidualBasedNewtonRaphsonStrategy */
} /* namespace Kratos. */

#endif /* KRATOS_UPWIND_RESIDUALBASED_NEWTON_RAPHSON_STRATEGY defined */
