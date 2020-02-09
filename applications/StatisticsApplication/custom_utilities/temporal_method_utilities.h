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

#if !defined(KRATOS_TEMPORAL_METHOD_UTILITIES_H_INCLUDED)
#define KRATOS_TEMPORAL_METHOD_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes

// Application includes
#include "custom_utilities/method_utilities.h"

namespace Kratos
{
///@addtogroup RANSApplication
///@{

///@name Kratos Globals
///@{

namespace TemporalMethodsUtilities
{
template <typename TContainerType, typename TContainerItemType, template <typename T> typename TDataRetrievalFunctor, template <typename T> typename TDataStorageFunctor, typename TDataType>
void InitializeVariables(TContainerType& rContainer,
                         const Variable<TDataType>& rOutputVariable,
                         const Variable<TDataType>& rReferenceVariable)
{
    if (rContainer.size() > 0)
    {
        const int number_of_items = rContainer.size();
#pragma omp parallel for
        for (int i = 0; i < number_of_items; ++i)
        {
            TContainerItemType& r_item = *(rContainer.begin() + i);
            const TDataType& r_reference_value =
                TDataRetrievalFunctor<TContainerItemType>()(r_item, rReferenceVariable);
            TDataType& r_output_value =
                TDataStorageFunctor<TContainerItemType>()(r_item, rOutputVariable);
            r_output_value = rOutputVariable.Zero();
            MethodsUtilities::DataTypeSizeInitializer<TDataType>(
                r_output_value, r_reference_value);
        }
    }
}

template <typename TContainerType, typename TContainerItemType, template <typename T> typename TDataStorageFunctor>
void InitializeVariables(TContainerType& rContainer,
                         const Variable<double>& rOutputVariable,
                         const double InitializerValue)
{
    if (rContainer.size() > 0)
    {
        const int number_of_items = rContainer.size();
#pragma omp parallel for
        for (int i = 0; i < number_of_items; ++i)
        {
            TContainerItemType& r_item = *(rContainer.begin() + i);
            double& r_output_value =
                TDataStorageFunctor<TContainerItemType>()(r_item, rOutputVariable);
            r_output_value = InitializerValue;
        }
    }
}

} // namespace TemporalMethodsUtilities
} // namespace Kratos
#endif // KRATOS_TEMPORAL_METHOD_UTILITIES_H_INCLUDED