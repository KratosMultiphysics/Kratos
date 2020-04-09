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
#include "factories/base_factory.h"

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
 * @tparam TStrategyType The strategy type
 * @tparam TCustomStrategyType The strategy type  (derived class)
 */
template <typename TStrategyType, typename TCustomStrategyType>
class StandardStrategyFactory
    : public BaseFactory<TStrategyType>
{
public:
    ///@name Type Definitions
    ///@{

    /// The definition of the strategy
    typedef TStrategyType StrategyType;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method returns the name of the class stored
     * @return The name of the class stored (as defined in settings)
     */
    std::string Name() const override
    {
        return TCustomStrategyType::Name();
    }

    /**
     * @brief This method is an auxiliar method to create a new strategy
     * @return The pointer to the strategy of interest
     */
    typename StrategyType::Pointer CreateClass(ModelPart& rModelPart, Kratos::Parameters Settings) const override
    {
        return Kratos::make_shared<TCustomStrategyType>(rModelPart, Settings);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "StandardStrategyFactory";
    }
};

///@}
///@name Input and output
///@{

/// output stream function
template <typename TStrategyType, typename TCustomStrategyType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const StandardStrategyFactory<TStrategyType, TCustomStrategyType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}
///@name Input and output

void RegisterStrategies();

///@}

}  // namespace Kratos.

#endif // KRATOS_STANDARD_STRATEGY_FACTORY_H_INCLUDED  defined
