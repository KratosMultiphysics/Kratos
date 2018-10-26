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

#if !defined(KRATOS_STRATEGY_FACTORY_H_INCLUDED )
#define  KRATOS_STRATEGY_FACTORY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/kratos_components.h"
#include "includes/shared_pointers.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "spaces/ublas_space.h"

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
 * @class StrategyFactory
 * @ingroup KratosCore
 * @brief Here we add the functions needed for the registration of strategies
 * @details Defines the base strategies factory
 * @author Vicente Mataix Ferrandiz
 * @tparam TSparseSpace The sparse space definition
 * @tparam TLocalSpace The dense space definition
 * @tparam TLinearSolver The linear solver type considered
 */
template< typename TSparseSpace, typename TLocalSpace, class TLinearSolver>
class StrategyFactory
{
public:
    ///@name Type Definitions
    ///@{

    /// The definition of the custom class
    typedef StrategyFactory<TSparseSpace,TLocalSpace, TLinearSolver> FactoryType;

    /// The definition of the strategies
    typedef SolvingStrategy<TSparseSpace,TLocalSpace, TLinearSolver> StrategyType;

    /// Pointer definition of StrategyFactory
    KRATOS_CLASS_POINTER_DEFINITION(StrategyFactory );

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method checks if the linear strategies is registered
     * @return True if registered, false otherwise
     */
    virtual bool Has(const std::string& rSolverType)
    {
        return KratosComponents< FactoryType >::Has( rSolverType );
    }

    /**
     * @brief This method creates a new strategies
     * @return The pointer to the strategies of interest
     */
    virtual typename StrategyType::Pointer Create(ModelPart& rModelPart, Kratos::Parameters Settings)
    {
        const std::string& strategy_type = Settings["strategy_type"].GetString();
        if(Has( strategy_type ) == false) {
            KRATOS_ERROR << "Trying to construct a strategies with type strategy_type= " << strategy_type << std::endl <<
                            "Which does not exist. The list of available options (for currently loaded applications) is: " << std::endl <<
                         KratosComponents< FactoryType >() << std::endl;
        }
        const auto& aux = KratosComponents< FactoryType >::Get( strategy_type );
        return aux.CreateStrategy(rModelPart, Settings);
    }
    ///@}
protected:
    ///@name Protected Operators
    ///@{

    /**
     * @brief This method is an auxiliar method to create a new strategies
     * @return The pointer to the strategies of interest
     */
    virtual typename StrategyType::Pointer CreateStrategy(ModelPart& rModelPart, Kratos::Parameters Settings)  const
    {
        KRATOS_ERROR << "Calling the base class StrategyFactory" << std::endl;
    }
};

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
template< typename TSparseSpace, typename TLocalSpace, class TLinearSolver>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const StrategyFactory<TSparseSpace, TLocalSpace, TLinearSolver>& rThis)
{
    rOStream << "StrategyFactory" << std::endl;

    return rOStream;
}

///@}

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSparseSpaceType;
typedef LinearSolver<SparseSpaceType,LocalSparseSpaceType> LinearSolverType;

typedef StrategyFactory<SparseSpaceType, LocalSparseSpaceType, LinearSolverType> StrategyFactoryType;

#ifdef KRATOS_REGISTER_STRATEGY
#undef KRATOS_REGISTER_STRATEGY
#endif
#define KRATOS_REGISTER_STRATEGY(name, reference) \
    KratosComponents<StrategyFactoryType>::Add(name, reference);

}  // namespace Kratos.

#endif // KRATOS_STRATEGY_FACTORY_H_INCLUDED  defined
