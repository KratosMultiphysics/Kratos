//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_SPATIAL_VARIANCE_METHOD_H_INCLUDED)
#define KRATOS_SPATIAL_VARIANCE_METHOD_H_INCLUDED

// System includes
#include <cmath>

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_utilities/method_utilities.h"

namespace Kratos
{
///@addtogroup RANSApplication
///@{

///@name Kratos Globals
///@{

namespace SpatialMethods
{
template <typename TContainerType, typename TContainerItem, template <typename T> typename TDataRetrievalFunctor>
class ContainerSpatialMethods
{
public:
    template <typename TDataType>
    void static CalculateSum(TDataType& rSum,
                             const TContainerType& rContainer,
                             const Variable<TDataType>& rVariable)
    {
        rSum = rVariable.Zero();
        MethodsUtilities::DataTypeSizeInitializer(
            rSum, TDataRetrievalFunctor<TContainerItem>()(*rContainer.begin(), rVariable));
        if (rContainer.size() > 0)
        {
            for (const auto& r_item : rContainer)
            {
                const TDataType& current_value =
                    TDataRetrievalFunctor<TContainerItem>()(r_item, rVariable);
                MethodsUtilities::DataTypeSizeChecker(current_value, rSum);
                rSum += current_value;
            }
        }
    }

    // special overloaded method for flags
    void static CalculateSum(int& rSum, const TContainerType& rContainer, const Flags& rVariable)
    {
        rSum = 0;
        if (rContainer.size() > 0)
        {
            for (const auto& r_item : rContainer)
            {
                rSum += r_item.Is(rVariable);
            }
        }
    }

    template <typename TDataType>
    void static CalculateMean(TDataType& rMean,
                              const TContainerType& rContainer,
                              const Variable<TDataType>& rVariable)
    {
        if (rContainer.size() > 0)
        {
            CalculateSum<TDataType>(rMean, rContainer, rVariable);
            rMean *= (1.0 / static_cast<double>(rContainer.size()));
        }
    }

    template <typename TDataType>
    void static CalculateVariance(TDataType& rMean,
                                  TDataType& rVariance,
                                  const TContainerType& rContainer,
                                  const Variable<TDataType>& rVariable)
    {
        CalculateMean<TContainerType, TContainerItem, TDataType, TDataRetrievalFunctor>(
            rMean, rContainer, rVariable);

        rVariance = rVariable.Zero();
        MethodsUtilities::DataTypeSizeInitializer(
            rVariance, TDataRetrievalFunctor<TContainerItem>()(*rContainer.begin(), rVariable));

        if (rContainer.size() > 0)
        {
            for (const auto& r_item : rContainer)
            {
                const TDataType& current_value =
                    TDataRetrievalFunctor<TContainerItem>()(r_item, rVariable);
                MethodsUtilities::DataTypeSizeChecker(current_value, rVariance);

                rVariance += MethodsUtilities::RaiseToPower(current_value, 2);
            }
            rVariance *= (1.0 / static_cast<double>(rContainer.size()));
            rVariance -= MethodsUtilities::RaiseToPower(rMean, 2);
        }
    }
};
} // namespace SpatialMethods

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

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_SPATIAL_VARIANCE_METHOD_H_INCLUDED defined
