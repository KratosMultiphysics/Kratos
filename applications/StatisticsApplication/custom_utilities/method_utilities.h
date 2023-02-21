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
#include <functional>

// External includes

// Project includes
#include "includes/model_part.h"

// Application includes

#ifndef ADD_TEMPORAL_VALUE_METHOD_ONE_OUTPUT_VARIABLE_OBJECT
#define ADD_TEMPORAL_VALUE_METHOD_ONE_OUTPUT_VARIABLE_OBJECT(                                \
    model_part, norm_type, input_variable, echo_level, output_variable, object_list, method) \
    {                                                                                        \
        if (KratosComponents<Variable<double>>::Has(input_variable))                         \
        {                                                                                    \
            const Variable<double>& r_variable_input =                                       \
                KratosComponents<Variable<double>>::Get(input_variable);                     \
            const Variable<double>& r_variable_output =                                      \
                KratosComponents<Variable<double>>::Get(output_variable);                    \
            object_list.push_back(std::make_shared<method<double>>(                          \
                model_part, norm_type, r_variable_input, echo_level, r_variable_output));    \
        }                                                                                    \
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(input_variable))       \
        {                                                                                    \
            const Variable<array_1d<double, 3>>& r_variable_input =                          \
                KratosComponents<Variable<array_1d<double, 3>>>::Get(input_variable);        \
            const Variable<array_1d<double, 3>>& r_variable_output =                         \
                KratosComponents<Variable<array_1d<double, 3>>>::Get(output_variable);       \
            object_list.push_back(std::make_shared<method<array_1d<double, 3>>>(             \
                model_part, norm_type, r_variable_input, echo_level, r_variable_output));    \
        }                                                                                    \
        else if (KratosComponents<Variable<Vector>>::Has(input_variable))                    \
        {                                                                                    \
            const Variable<Vector>& r_variable_input =                                       \
                KratosComponents<Variable<Vector>>::Get(input_variable);                     \
            const Variable<Vector>& r_variable_output =                                      \
                KratosComponents<Variable<Vector>>::Get(output_variable);                    \
            object_list.push_back(std::make_shared<method<Vector>>(                          \
                model_part, norm_type, r_variable_input, echo_level, r_variable_output));    \
        }                                                                                    \
        else if (KratosComponents<Variable<Matrix>>::Has(input_variable))                    \
        {                                                                                    \
            const Variable<Matrix>& r_variable_input =                                       \
                KratosComponents<Variable<Matrix>>::Get(input_variable);                     \
            const Variable<Matrix>& r_variable_output =                                      \
                KratosComponents<Variable<Matrix>>::Get(output_variable);                    \
            object_list.push_back(std::make_shared<method<Matrix>>(                          \
                model_part, norm_type, r_variable_input, echo_level, r_variable_output));    \
        }                                                                                    \
        else                                                                                 \
        {                                                                                    \
            KRATOS_ERROR                                                                     \
                << "Input variable not found in Double, Array3D, Vector or "                 \
                   "Matrix variables list. [ input_variable = "                              \
                << input_variable << " ]\n";                                                 \
        }                                                                                    \
    }
#endif

#ifndef ADD_TEMPORAL_NORM_METHOD_ONE_OUTPUT_VARIABLE_OBJECT
#define ADD_TEMPORAL_NORM_METHOD_ONE_OUTPUT_VARIABLE_OBJECT(                                 \
    model_part, norm_type, input_variable, echo_level, output_variable, object_list, method) \
    {                                                                                        \
        if (KratosComponents<Variable<double>>::Has(input_variable))                         \
        {                                                                                    \
            const Variable<double>& r_variable_input =                                       \
                KratosComponents<Variable<double>>::Get(input_variable);                     \
            const Variable<double>& r_variable_output =                                      \
                KratosComponents<Variable<double>>::Get(output_variable);                    \
            object_list.push_back(std::make_shared<method<double>>(                          \
                model_part, norm_type, r_variable_input, echo_level, r_variable_output));    \
        }                                                                                    \
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(input_variable))       \
        {                                                                                    \
            const Variable<array_1d<double, 3>>& r_variable_input =                          \
                KratosComponents<Variable<array_1d<double, 3>>>::Get(input_variable);        \
            const Variable<double>& r_variable_output =                                      \
                KratosComponents<Variable<double>>::Get(output_variable);                    \
            object_list.push_back(std::make_shared<method<array_1d<double, 3>>>(             \
                model_part, norm_type, r_variable_input, echo_level, r_variable_output));    \
        }                                                                                    \
        else if (KratosComponents<Variable<Vector>>::Has(input_variable))                    \
        {                                                                                    \
            const Variable<Vector>& r_variable_input =                                       \
                KratosComponents<Variable<Vector>>::Get(input_variable);                     \
            const Variable<double>& r_variable_output =                                      \
                KratosComponents<Variable<double>>::Get(output_variable);                    \
            object_list.push_back(std::make_shared<method<Vector>>(                          \
                model_part, norm_type, r_variable_input, echo_level, r_variable_output));    \
        }                                                                                    \
        else if (KratosComponents<Variable<Matrix>>::Has(input_variable))                    \
        {                                                                                    \
            const Variable<Matrix>& r_variable_input =                                       \
                KratosComponents<Variable<Matrix>>::Get(input_variable);                     \
            const Variable<double>& r_variable_output =                                      \
                KratosComponents<Variable<double>>::Get(output_variable);                    \
            object_list.push_back(std::make_shared<method<Matrix>>(                          \
                model_part, norm_type, r_variable_input, echo_level, r_variable_output));    \
        }                                                                                    \
        else                                                                                 \
        {                                                                                    \
            KRATOS_ERROR                                                                     \
                << "Input variable not found in Double, Array3D, Vector or "                 \
                   "Matrix variables list. [ input_variable = "                              \
                << input_variable << " ]\n";                                                 \
        }                                                                                    \
    }
#endif

#ifndef ADD_TEMPORAL_VALUE_METHOD_TWO_OUTPUT_VARIABLE_OBJECT
#define ADD_TEMPORAL_VALUE_METHOD_TWO_OUTPUT_VARIABLE_OBJECT(                            \
    model_part, norm_type, input_variable, echo_level, output_variable_1,                \
    output_variable_2, object_list, method)                                              \
    {                                                                                    \
        if (KratosComponents<Variable<double>>::Has(input_variable))                     \
        {                                                                                \
            const Variable<double>& r_variable_input =                                   \
                KratosComponents<Variable<double>>::Get(input_variable);                 \
            const Variable<double>& r_variable_output_1 =                                \
                KratosComponents<Variable<double>>::Get(output_variable_1);              \
            const Variable<double>& r_variable_output_2 =                                \
                KratosComponents<Variable<double>>::Get(output_variable_2);              \
            object_list.push_back(std::make_shared<method<double>>(                      \
                model_part, norm_type, r_variable_input, echo_level,                     \
                r_variable_output_1, r_variable_output_2));                              \
        }                                                                                \
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(input_variable))   \
        {                                                                                \
            const Variable<array_1d<double, 3>>& r_variable_input =                      \
                KratosComponents<Variable<array_1d<double, 3>>>::Get(input_variable);    \
            const Variable<array_1d<double, 3>>& r_variable_output_1 =                   \
                KratosComponents<Variable<array_1d<double, 3>>>::Get(output_variable_1); \
            const Variable<array_1d<double, 3>>& r_variable_output_2 =                   \
                KratosComponents<Variable<array_1d<double, 3>>>::Get(output_variable_2); \
            object_list.push_back(std::make_shared<method<array_1d<double, 3>>>(         \
                model_part, norm_type, r_variable_input, echo_level,                     \
                r_variable_output_1, r_variable_output_2));                              \
        }                                                                                \
        else if (KratosComponents<Variable<Vector>>::Has(input_variable))                \
        {                                                                                \
            const Variable<Vector>& r_variable_input =                                   \
                KratosComponents<Variable<Vector>>::Get(input_variable);                 \
            const Variable<Vector>& r_variable_output_1 =                                \
                KratosComponents<Variable<Vector>>::Get(output_variable_1);              \
            const Variable<Vector>& r_variable_output_2 =                                \
                KratosComponents<Variable<Vector>>::Get(output_variable_2);              \
            object_list.push_back(std::make_shared<method<Vector>>(                      \
                model_part, norm_type, r_variable_input, echo_level,                     \
                r_variable_output_1, r_variable_output_2));                              \
        }                                                                                \
        else if (KratosComponents<Variable<Matrix>>::Has(input_variable))                \
        {                                                                                \
            const Variable<Matrix>& r_variable_input =                                   \
                KratosComponents<Variable<Matrix>>::Get(input_variable);                 \
            const Variable<Matrix>& r_variable_output_1 =                                \
                KratosComponents<Variable<Matrix>>::Get(output_variable_1);              \
            const Variable<Matrix>& r_variable_output_2 =                                \
                KratosComponents<Variable<Matrix>>::Get(output_variable_2);              \
            object_list.push_back(std::make_shared<method<Matrix>>(                      \
                model_part, norm_type, r_variable_input, echo_level,                     \
                r_variable_output_1, r_variable_output_2));                              \
        }                                                                                \
        else                                                                             \
        {                                                                                \
            KRATOS_ERROR                                                                 \
                << "Input variable not found in Double, Array3D, Vector or "             \
                   "Matrix variables list. [ input_variable = "                          \
                << input_variable << " ]\n";                                             \
        }                                                                                \
    }
#endif

#ifndef ADD_TEMPORAL_NORM_METHOD_TWO_OUTPUT_VARIABLE_OBJECT
#define ADD_TEMPORAL_NORM_METHOD_TWO_OUTPUT_VARIABLE_OBJECT(                           \
    model_part, norm_type, input_variable, echo_level, output_variable_1,              \
    output_variable_2, object_list, method)                                            \
    {                                                                                  \
        if (KratosComponents<Variable<double>>::Has(input_variable))                   \
        {                                                                              \
            const Variable<double>& r_variable_input =                                 \
                KratosComponents<Variable<double>>::Get(input_variable);               \
            const Variable<double>& r_variable_output_1 =                              \
                KratosComponents<Variable<double>>::Get(output_variable_1);            \
            const Variable<double>& r_variable_output_2 =                              \
                KratosComponents<Variable<double>>::Get(output_variable_2);            \
            object_list.push_back(std::make_shared<method<double>>(                    \
                model_part, norm_type, r_variable_input, echo_level,                   \
                r_variable_output_1, r_variable_output_2));                            \
        }                                                                              \
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(input_variable)) \
        {                                                                              \
            const Variable<array_1d<double, 3>>& r_variable_input =                    \
                KratosComponents<Variable<array_1d<double, 3>>>::Get(input_variable);  \
            const Variable<double>& r_variable_output_1 =                              \
                KratosComponents<Variable<double>>::Get(output_variable_1);            \
            const Variable<double>& r_variable_output_2 =                              \
                KratosComponents<Variable<double>>::Get(output_variable_2);            \
            object_list.push_back(std::make_shared<method<array_1d<double, 3>>>(       \
                model_part, norm_type, r_variable_input, echo_level,                   \
                r_variable_output_1, r_variable_output_2));                            \
        }                                                                              \
        else if (KratosComponents<Variable<Vector>>::Has(input_variable))              \
        {                                                                              \
            const Variable<Vector>& r_variable_input =                                 \
                KratosComponents<Variable<Vector>>::Get(input_variable);               \
            const Variable<double>& r_variable_output_1 =                              \
                KratosComponents<Variable<double>>::Get(output_variable_1);            \
            const Variable<double>& r_variable_output_2 =                              \
                KratosComponents<Variable<double>>::Get(output_variable_2);            \
            object_list.push_back(std::make_shared<method<Vector>>(                    \
                model_part, norm_type, r_variable_input, echo_level,                   \
                r_variable_output_1, r_variable_output_2));                            \
        }                                                                              \
        else if (KratosComponents<Variable<Matrix>>::Has(input_variable))              \
        {                                                                              \
            const Variable<Matrix>& r_variable_input =                                 \
                KratosComponents<Variable<Matrix>>::Get(input_variable);               \
            const Variable<double>& r_variable_output_1 =                              \
                KratosComponents<Variable<double>>::Get(output_variable_1);            \
            const Variable<double>& r_variable_output_2 =                              \
                KratosComponents<Variable<double>>::Get(output_variable_2);            \
            object_list.push_back(std::make_shared<method<Matrix>>(                    \
                model_part, norm_type, r_variable_input, echo_level,                   \
                r_variable_output_1, r_variable_output_2));                            \
        }                                                                              \
        else                                                                           \
        {                                                                              \
            KRATOS_ERROR                                                               \
                << "Input variable not found in Double, Array3D, Vector or "           \
                   "Matrix variables list. [ input_variable = "                        \
                << input_variable << " ]\n";                                           \
        }                                                                              \
    }
#endif

namespace Kratos
{
///@addtogroup RANSApplication
///@{

///@name Kratos Globals
///@{

namespace MethodUtilities
{
using NodeType = ModelPart::NodeType;
using ConditionType = ModelPart::ConditionType;
using ElementType = ModelPart::ElementType;

using NodesContainerType = ModelPart::NodesContainerType;
using ConditionsContainerType = ModelPart::ConditionsContainerType;
using ElementsContainerType = ModelPart::ElementsContainerType;

template <class TDataType>
KRATOS_API(STATISTICS_APPLICATION)
TDataType RaiseToPower(const TDataType& rData, const double Power);

template <class TContainerItemType>
class NonHistoricalDataValueRetrievalFunctor
{
public:
    template <class TDataType>
    TDataType& operator()(TContainerItemType& rDataItem, const Variable<TDataType>& rVariable) const
    {
        return rDataItem.GetValue(rVariable);
    }

    template <class TDataType>
    TDataType operator()(const TContainerItemType& rDataItem, const Variable<TDataType>& rVariable) const
    {
        return rDataItem.GetValue(rVariable);
    }

    template <class TDataType>
    void operator()(TContainerItemType& rDataItem, const Variable<TDataType>& rVariable, const TDataType& rValue) const
    {
        rDataItem.SetValue(rVariable, rValue);
    }
};

template <class TContainerItemType>
class HistoricalDataValueRetrievalFunctor
{
public:
    template <class TDataType>
    TDataType& operator()(TContainerItemType& rDataItem, const Variable<TDataType>& rVariable) const
    {
        KRATOS_TRY

        return rDataItem.FastGetSolutionStepValue(rVariable);

        KRATOS_CATCH("");
    }

    template <class TDataType>
    TDataType operator()(const TContainerItemType& rDataItem, const Variable<TDataType>& rVariable) const
    {
        KRATOS_TRY

        return rDataItem.FastGetSolutionStepValue(rVariable);

        KRATOS_CATCH("");
    }

    template <class TDataType>
    void operator()(TContainerItemType& rDataItem, const Variable<TDataType>& rVariable, const TDataType& rValue) const
    {
        KRATOS_TRY

        rDataItem.FastGetSolutionStepValue(rVariable) = rValue;

        KRATOS_CATCH("");
    }
};

KRATOS_API(STATISTICS_APPLICATION)
double GetDoubleValue(const std::string& rInput);

KRATOS_API(STATISTICS_APPLICATION)
int GetIntegerValue(const std::string& rInput);

KRATOS_API(STATISTICS_APPLICATION)
void SplitString(std::string& rOutput1, std::string& rOutput2, const std::string& rInput);

template <class TDataType>
KRATOS_API(STATISTICS_APPLICATION)
void DataTypeSizeInitializer(TDataType& rData, const TDataType& rReferenceData);

template <class TDataType>
KRATOS_API(STATISTICS_APPLICATION)
void DataTypeSizeChecker(const TDataType& rData, const TDataType& rReferenceData);

template <class TContainerType>
KRATOS_API(STATISTICS_APPLICATION)
TContainerType& GetLocalDataContainer(ModelPart& rModelPart);

template <class TContainerType>
KRATOS_API(STATISTICS_APPLICATION)
const TContainerType& GetLocalDataContainer(const ModelPart& rModelPart);

template <class TContainerType>
KRATOS_API(STATISTICS_APPLICATION)
TContainerType& GetDataContainer(ModelPart& rModelPart);

template <class TContainerType>
KRATOS_API(STATISTICS_APPLICATION)
const TContainerType& GetDataContainer(const ModelPart& rModelPart);

template <class TDataType>
KRATOS_API(STATISTICS_APPLICATION)
const std::function<double(const TDataType&)> GetNormMethod(
    const Variable<TDataType>& rVariable, const std::string& rNormType);

template <class TDataType>
KRATOS_API(STATISTICS_APPLICATION)
std::string GetVariableTypeName();

template <class TDataType>
KRATOS_API(STATISTICS_APPLICATION)
void CheckVariableType(const std::vector<std::string>& rVariableNamesList);

KRATOS_API(STATISTICS_APPLICATION)
void CheckInputOutputVariables(
    const std::vector<std::string>& rInputVariableNamesList,
    const std::vector<std::string>& rOutputVariableNamesList);

KRATOS_API(STATISTICS_APPLICATION)
std::vector<double> SortSortedValuesList(const std::vector<std::vector<double>>& rValues);

} // namespace MethodUtilities

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
