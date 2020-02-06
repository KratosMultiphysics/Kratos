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

#if !defined(KRATOS_METHOD_UTILITIES_H_INCLUDED)
#define KRATOS_METHOD_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"

// Application includes

namespace Kratos
{
///@addtogroup RANSApplication
///@{

///@name Kratos Globals
///@{

namespace MethodsUtilities
{
using NodeType = ModelPart::NodeType;
using ConditionType = ModelPart::ConditionType;
using ElementType = ModelPart::ElementType;

using NodesContainerType = ModelPart::NodesContainerType;
using ConditionsContainerType = ModelPart::ConditionsContainerType;
using ElementsContainerType = ModelPart::ElementsContainerType;

template <typename TDataType>
TDataType RaiseToPower(const TDataType& rData, const double Power);

template <typename TContainerItemType>
class NonHistoricalDataValueRetrievalFunctor
{
public:
    template <typename TDataType>
    TDataType& operator()(TContainerItemType& rDataItem, const Variable<TDataType>& rVariable)
    {
        return rDataItem.GetValue(rVariable);
    }

    template <typename TDataType>
    TDataType operator()(const TContainerItemType& rDataItem,
                         const Variable<TDataType>& rVariable)
    {
        return rDataItem.GetValue(rVariable);
    }
};

template <typename TContainerItemType>
class HistoricalDataValueRetrievalFunctor
{
public:
    template <typename TDataType>
    TDataType& operator()(TContainerItemType& rDataItem, const Variable<TDataType>& rVariable)
    {
        KRATOS_TRY

        return rDataItem.FastGetSolutionStepValue(rVariable);

        KRATOS_CATCH("");
    }

    template <typename TDataType>
    TDataType operator()(const TContainerItemType& rDataItem,
                         const Variable<TDataType>& rVariable)
    {
        KRATOS_TRY

        return rDataItem.FastGetSolutionStepValue(rVariable);

        KRATOS_CATCH("");
    }
};

template <typename TDataType>
void DataTypeSizeInitializer(TDataType& rData, const TDataType& rReferenceData);

template <typename TDataType>
void DataTypeSizeChecker(const TDataType& rData, const TDataType& rReferenceData);

template <typename TContainerType>
TContainerType& GetDataContainer(ModelPart& rModelPart);

template <typename TContainerType>
const TContainerType& GetDataContainer(const ModelPart& rModelPart);

} // namespace MethodsUtilities

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

#endif // KRATOS_METHOD_UTILITIES_H_INCLUDED defined
