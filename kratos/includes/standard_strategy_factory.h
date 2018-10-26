//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

#if !defined(KRATOS_STANDARD_STRATEGY_FACTORY_H_INCLUDED )
#define  KRATOS_STANDARD_STRATEGY_FACTORY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/strategy_factory.h"
#include "solving_strategies/strategies/solving_strategy.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class StandardStrategyFactory
 * @ingroup KratosCore
 * @brief Here we add the functions needed for the registration of strategys
 * @details Defines the standard strategy factory
 * @author Vicente Mataix Ferrandiz
 * @tparam TSparseSpace The sparse space definition
 * @tparam TLocalSpace The dense space definition
 * @tparam TLinearSolver The linear solver type considered
 * @tparam TStrategyType The strategy type
 */
template <typename TSparseSpace, typename TLocalSpace, class TLinearSolver, typename TStrategyType>
class StandardStrategyFactory
    : public StrategyFactory<TSparseSpace,TLocalSpace, TLinearSolver>
{
    ///@name Type Definitions
    ///@{

    /// The definition of the strategy
    typedef SolvingStrategy<TSparseSpace,TLocalSpace, TLinearSolver> StrategyType;

    ///@}
protected:

    ///@name Operations
    ///@{

    /**
     * @brief This method is an auxiliar method to create a new strategy
     * @return The pointer to the strategy of interest
     */
    typename StrategyType::Pointer CreateStrategy(ModelPart& rModelPart, Kratos::Parameters Settings) const override
    {
        return typename StrategyType::Pointer(new TStrategyType(rModelPart, Settings));
    }

    ///@}
};

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
template <typename TSparseSpace, typename TLocalSpace, class TLinearSolver, typename TStrategyType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const StandardStrategyFactory<TSparseSpace,TLocalSpace, TLinearSolver, TStrategyType>& rThis)
{
    rOStream << "StandardStrategyFactory" << std::endl;
    return rOStream;
}

///@}
///@name Input and output

void RegisterStrategies();

///@}

}  // namespace Kratos.

#endif // KRATOS_STANDARD_STRATEGY_FACTORY_H_INCLUDED  defined
