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

#if !defined(KRATOS_TEMPORAL_SUM_METHOD_H_INCLUDED)
#define KRATOS_TEMPORAL_SUM_METHOD_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes
#include "custom_methods/temporal_method.h"
#include "custom_utilities/method_utilities.h"
#include "custom_utilities/temporal_method_utilities.h"

namespace Kratos
{
///@addtogroup StatisticsApplication
///@{

///@name Kratos Globals
///@{

namespace TemporalMethods
{
template <typename TContainerType, typename TContainerItemType, template <typename T> typename TDataRetrievalFunctor, template <typename T> typename TDataStorageFunctor>
class TemporalSumMethod : public TemporalMethod
{
public:
    using BaseType = TemporalMethod;

    KRATOS_CLASS_POINTER_DEFINITION(TemporalSumMethod);

    TemporalSumMethod(ModelPart& rModelPart) : BaseType(rModelPart)
    {
    }

    // value method
    template <typename TDataType>
    void CalculateValueStatistics(const Variable<TDataType>& rOutputVariable,
                                  const Variable<TDataType>& rInputVariable,
                                  const double DeltaTime)
    {
        TContainerType& r_container =
            MethodsUtilities::GetDataContainer<TContainerType>(this->GetModelPart());

        const int number_of_items = r_container.size();
#pragma omp parallel for
        for (int i = 0; i < number_of_items; ++i)
        {
            TContainerItemType& r_item = *(r_container.begin() + i);
            const TDataType& r_input_value =
                TDataRetrievalFunctor<TContainerItemType>()(r_item, rInputVariable);
            TDataType& r_output_value =
                TDataStorageFunctor<TContainerItemType>()(r_item, rOutputVariable);
            MethodsUtilities::DataTypeSizeChecker(r_input_value, r_output_value);

            CalculateSum<TDataType>(r_output_value, r_input_value, DeltaTime,
                                    this->GetTotalTime());
        }
    }

    // Norm method
    template <typename TDataType>
    void CalculateNormStatistics(const std::string& rNormType,
                                 const Variable<double>& rOutputVariable,
                                 const Variable<TDataType>& rInputVariable,
                                 const double DeltaTime)
    {
        TContainerType& r_container =
            MethodsUtilities::GetDataContainer<TContainerType>(this->GetModelPart());

        const auto& norm_method = MethodsUtilities::GetNormMethod(rInputVariable, rNormType);

        const int number_of_items = r_container.size();
#pragma omp parallel for
        for (int i = 0; i < number_of_items; ++i)
        {
            TContainerItemType& r_item = *(r_container.begin() + i);
            const TDataType& r_input_value =
                TDataRetrievalFunctor<TContainerItemType>()(r_item, rInputVariable);
            const double input_norm_value = norm_method(r_input_value);
            double& r_output_value =
                TDataStorageFunctor<TContainerItemType>()(r_item, rOutputVariable);

            CalculateSum<double>(r_output_value, input_norm_value, DeltaTime,
                                 this->GetTotalTime());
        }
    }

    // value output variable initialization
    template <typename TDataType>
    void InitializeValueStatisticsVariables(const Variable<TDataType>& rOutputVariable,
                                            const Variable<TDataType>& rInputVariable)
    {
        TContainerType& r_container =
            MethodsUtilities::GetDataContainer<TContainerType>(this->GetModelPart());

        auto& initializer_method =
            TemporalMethodsUtilities::InitializeVariables<TContainerType, TContainerItemType, TDataRetrievalFunctor,
                                                          TDataStorageFunctor, TDataType>;
        initializer_method(r_container, rOutputVariable, rInputVariable);
    }

    // norm output variable initialization
    void InitializeNormStatisticsVariables(const Variable<double>& rOutputVariable)
    {
        TContainerType& r_container =
            MethodsUtilities::GetDataContainer<TContainerType>(this->GetModelPart());

        auto& initializer_method =
            TemporalMethodsUtilities::InitializeVariables<TContainerType, TContainerItemType, TDataStorageFunctor>;
        initializer_method(r_container, rOutputVariable, 0.0);
    }

private:
    template <typename TDataType>
    void CalculateSum(TDataType& rSum,
                      const TDataType& rNewDataPoint,
                      const double DeltaTime,
                      const double TotalTime)
    {
        rSum = (rSum + rNewDataPoint);
    }
};
} // namespace TemporalMethods
} // namespace Kratos

#endif // KRATOS_TEMPORAL_SUM_METHOD_H_INCLUDED