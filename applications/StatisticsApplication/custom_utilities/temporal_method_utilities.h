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

namespace TemporalMethodUtilities
{
template <class TContainerType, class TContainerItemType, template <class T> class TDataRetrievalFunctor, template <class T> class TDataStorageFunctor, class TDataType>
void InitializeVariables(
    TContainerType& rContainer,
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
            TDataType output_value = rOutputVariable.Zero();
            MethodUtilities::DataTypeSizeInitializer<TDataType>(output_value, r_reference_value);
            TDataStorageFunctor<TContainerItemType>()(r_item, rOutputVariable, output_value);
        }
    }
}

template <class TContainerType, class TContainerItemType, template <class T> class TDataStorageFunctor>
void InitializeVariables(TContainerType& rContainer, const Variable<double>& rOutputVariable, const double InitializerValue)
{
    if (rContainer.size() > 0)
    {
        const int number_of_items = rContainer.size();
#pragma omp parallel for
        for (int i = 0; i < number_of_items; ++i)
        {
            TContainerItemType& r_item = *(rContainer.begin() + i);
            TDataStorageFunctor<TContainerItemType>()(r_item, rOutputVariable, InitializerValue);
        }
    }
}

} // namespace TemporalMethodUtilities
} // namespace Kratos
#endif // KRATOS_TEMPORAL_METHOD_UTILITIES_H_INCLUDED