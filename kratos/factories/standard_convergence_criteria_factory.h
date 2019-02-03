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

#if !defined(KRATOS_STANDARD_CONVERGENCE_CRITERIA_FACTORY_H_INCLUDED )
#define  KRATOS_STANDARD_CONVERGENCE_CRITERIA_FACTORY_H_INCLUDED

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
 * @class StandardConvergenceCriteriaFactory
 * @ingroup KratosCore
 * @brief Here we add the functions needed for the registration of convergence criterias
 * @details Defines the standard convergence criteria factory
 * @author Vicente Mataix Ferrandiz
 * @tparam TConvergenceCriteriaType The convergence criteria type
 * @tparam TCustomConvergenceCriteriaType The convergence criteria type (derived class)
 */
template <typename TConvergenceCriteriaType, typename TCustomConvergenceCriteriaType>
class StandardConvergenceCriteriaFactory
    : public BaseFactory<TConvergenceCriteriaType>
{
    ///@name Type Definitions
    ///@{

    /// The definition of the convergence criteria
    typedef TConvergenceCriteriaType ConvergenceCriteriaType;

    ///@}
protected:

    ///@name Operations
    ///@{

    /**
     * @brief This method is an auxiliar method to create a new convergence criteria
     * @return The pointer to the convergence criteria of interest
     */
    typename ConvergenceCriteriaType::Pointer CreateClass(Kratos::Parameters Settings) const override
    {
        return typename ConvergenceCriteriaType::Pointer(new TCustomConvergenceCriteriaType(Settings));
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
template <typename TConvergenceCriteriaType, typename TCustomConvergenceCriteriaType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const StandardConvergenceCriteriaFactory<TConvergenceCriteriaType, TCustomConvergenceCriteriaType>& rThis)
{
    rOStream << "StandardConvergenceCriteriaFactory" << std::endl;
    return rOStream;
}

///@}
///@name Input and output

void RegisterConvergenceCriterias();

///@}

}  // namespace Kratos.

#endif // KRATOS_STANDARD_CONVERGENCE_CRITERIA_FACTORY_H_INCLUDED  defined
