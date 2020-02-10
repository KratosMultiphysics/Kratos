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

#if !defined(KRATOS_TEMPORAL_VARIANCE_METHOD_H_INCLUDED)
#define KRATOS_TEMPORAL_VARIANCE_METHOD_H_INCLUDED

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
class TemporalVarianceMethod
{
public:
    template <typename TDataType>
    class ValueMethod : public TemporalMethod
    {
    public:
        KRATOS_CLASS_POINTER_DEFINITION(ValueMethod);

        ValueMethod(ModelPart& rModelPart,
                    const std::string& rNormType,
                    const Variable<TDataType>& rInputVariable,
                    const Variable<TDataType>& rOutputMeanVariable,
                    const Variable<TDataType>& rOutputVarianceVariable)
            : TemporalMethod(rModelPart),
              mrInputVariable(rInputVariable),
              mrOutputMeanVariable(rOutputMeanVariable),
              mrOutputVarianceVariable(rOutputVarianceVariable)
        {
        }

        void CalculateStatistics(const double DeltaTime) override
        {
            TContainerType& r_container =
                MethodsUtilities::GetDataContainer<TContainerType>(this->GetModelPart());

            const int number_of_items = r_container.size();
#pragma omp parallel for
            for (int i = 0; i < number_of_items; ++i)
            {
                TContainerItemType& r_item = *(r_container.begin() + i);
                const TDataType& r_input_value =
                    TDataRetrievalFunctor<TContainerItemType>()(r_item, mrInputVariable);
                TDataType& r_output_mean_value =
                    TDataStorageFunctor<TContainerItemType>()(r_item, mrOutputMeanVariable);
                TDataType& r_output_variance_value =
                    TDataStorageFunctor<TContainerItemType>()(r_item, mrOutputVarianceVariable);

                MethodsUtilities::DataTypeSizeChecker(r_input_value, r_output_mean_value);
                MethodsUtilities::DataTypeSizeChecker(r_input_value, r_output_variance_value);

                TemporalVarianceMethod::CalculateMeanAndVariance<TDataType>(
                    r_output_mean_value, r_output_variance_value, r_input_value,
                    DeltaTime, this->GetTotalTime());
            }

            TemporalMethod::CalculateStatistics(DeltaTime);
        }

        void InitializeStatisticsVariables() override
        {
            TContainerType& r_container =
                MethodsUtilities::GetDataContainer<TContainerType>(this->GetModelPart());

            auto& initializer_method =
                TemporalMethodsUtilities::InitializeVariables<TContainerType, TContainerItemType, TDataRetrievalFunctor,
                                                              TDataStorageFunctor, TDataType>;
            initializer_method(r_container, mrOutputMeanVariable, mrInputVariable);
            initializer_method(r_container, mrOutputVarianceVariable, mrInputVariable);
        }

    private:
        const Variable<TDataType>& mrInputVariable;
        const Variable<TDataType>& mrOutputMeanVariable;
        const Variable<TDataType>& mrOutputVarianceVariable;
    };

    template <typename TDataType>
    class NormMethod : public TemporalMethod
    {
    public:
        KRATOS_CLASS_POINTER_DEFINITION(NormMethod);

        NormMethod(ModelPart& rModelPart,
                   const std::string& rNormType,
                   const Variable<TDataType>& rInputVariable,
                   const Variable<double>& rOutputMeanVariable,
                   const Variable<double>& rOutputVarianceVariable)
            : TemporalMethod(rModelPart),
              mNormType(rNormType),
              mrInputVariable(rInputVariable),
              mrOutputMeanVariable(rOutputMeanVariable),
              mrOutputVarianceVariable(rOutputVarianceVariable)
        {
        }

        void CalculateStatistics(const double DeltaTime) override
        {
            TContainerType& r_container =
                MethodsUtilities::GetDataContainer<TContainerType>(this->GetModelPart());

            const auto& norm_method =
                MethodsUtilities::GetNormMethod(mrInputVariable, mNormType);

            const int number_of_items = r_container.size();
#pragma omp parallel for
            for (int i = 0; i < number_of_items; ++i)
            {
                TContainerItemType& r_item = *(r_container.begin() + i);
                const TDataType& r_input_value =
                    TDataRetrievalFunctor<TContainerItemType>()(r_item, mrInputVariable);
                const double input_norm_value = norm_method(r_input_value);
                double& r_output_mean_value =
                    TDataStorageFunctor<TContainerItemType>()(r_item, mrOutputMeanVariable);
                double& r_output_variance_value =
                    TDataStorageFunctor<TContainerItemType>()(r_item, mrOutputVarianceVariable);

                TemporalVarianceMethod::CalculateMeanAndVariance<double>(
                    r_output_mean_value, r_output_variance_value,
                    input_norm_value, DeltaTime, this->GetTotalTime());
            }

            TemporalMethod::CalculateStatistics(DeltaTime);
        }

        // norm output variable initialization
        void InitializeStatisticsVariables() override
        {
            TContainerType& r_container =
                MethodsUtilities::GetDataContainer<TContainerType>(this->GetModelPart());

            auto& initializer_method =
                TemporalMethodsUtilities::InitializeVariables<TContainerType, TContainerItemType, TDataStorageFunctor>;
            initializer_method(r_container, mrOutputMeanVariable, 0.0);
            initializer_method(r_container, mrOutputVarianceVariable, 0.0);
        }

    private:
        const std::string mNormType;
        const Variable<TDataType>& mrInputVariable;
        const Variable<double>& mrOutputMeanVariable;
        const Variable<double>& mrOutputVarianceVariable;
    };

private:
    template <typename TDataType>
    void static CalculateMeanAndVariance(TDataType& rMean,
                                         TDataType& rVariance,
                                         const TDataType& rNewDataPoint,
                                         const double DeltaTime,
                                         const double TotalTime)
    {
        const double new_total_time = TotalTime + DeltaTime;
        const TDataType new_mean =
            (rMean * TotalTime + rNewDataPoint) * (1.0 / new_total_time);
        rVariance =
            ((rVariance + MethodsUtilities::RaiseToPower<TDataType>(new_mean - rMean, 2)) * TotalTime +
             MethodsUtilities::RaiseToPower<TDataType>(rNewDataPoint - new_mean, 2)) *
            (1.0 / new_total_time);
        rMean = new_mean;
    }
};
} // namespace TemporalMethods
} // namespace Kratos

#endif // KRATOS_TEMPORAL_VARIANCE_METHOD_H_INCLUDED