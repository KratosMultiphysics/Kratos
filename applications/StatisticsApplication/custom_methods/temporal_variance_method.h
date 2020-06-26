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
template <class TContainerType, class TContainerItemType, template <class T> class TDataRetrievalFunctor, template <class T> class TDataStorageFunctor>
class TemporalVarianceMethod
{
public:
    template <class TDataType>
    class ValueMethod : public TemporalMethod
    {
    public:
        KRATOS_CLASS_POINTER_DEFINITION(ValueMethod);

        ValueMethod(
            ModelPart& rModelPart,
            const std::string& rNormType,
            const Variable<TDataType>& rInputVariable,
            const int EchoLevel,
            const Variable<TDataType>& rOutputMeanVariable,
            const Variable<TDataType>& rOutputVarianceVariable)
            : TemporalMethod(rModelPart, EchoLevel),
              mrInputVariable(rInputVariable),
              mrOutputMeanVariable(rOutputMeanVariable),
              mrOutputVarianceVariable(rOutputVarianceVariable)
        {
            KRATOS_TRY

            KRATOS_ERROR_IF(rOutputMeanVariable == rOutputVarianceVariable) << "Same variable is given for mean and variance in value variance method with input variable "
                                                                            << rInputVariable
                                                                                   .Name()
                                                                            << ". Please provide two different variables. [ variable = "
                                                                            << rOutputMeanVariable
                                                                                   .Name()
                                                                            << " ].\n";

            KRATOS_CATCH("");
        }

        void CalculateStatistics() override
        {
            TContainerType& r_container =
                MethodUtilities::GetDataContainer<TContainerType>(this->GetModelPart());

            const double delta_time = this->GetDeltaTime();
            const double old_total_time = this->GetTotalTime();
            const double total_time = old_total_time + delta_time;

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

                MethodUtilities::DataTypeSizeChecker(r_input_value, r_output_mean_value);
                MethodUtilities::DataTypeSizeChecker(r_input_value, r_output_variance_value);

                TemporalVarianceMethod::CalculateMeanAndVariance<TDataType>(
                    r_output_mean_value, r_output_variance_value, r_input_value,
                    delta_time, old_total_time, total_time);
            }

            KRATOS_INFO_IF("TemporalValueVarianceMethod", this->GetEchoLevel() > 1)
                << "Calculated temporal value variance for "
                << mrInputVariable.Name() << " input variable with "
                << mrOutputMeanVariable.Name() << " mean variable and "
                << mrOutputVarianceVariable.Name() << " variance variable for "
                << this->GetModelPart().Name() << ".\n";
        }

        void InitializeStatisticsVariables() override
        {
            TContainerType& r_container =
                MethodUtilities::GetDataContainer<TContainerType>(this->GetModelPart());

            auto& initializer_method =
                TemporalMethodUtilities::InitializeVariables<TContainerType, TContainerItemType, TDataRetrievalFunctor, TDataStorageFunctor, TDataType>;
            initializer_method(r_container, mrOutputMeanVariable, mrInputVariable);
            initializer_method(r_container, mrOutputVarianceVariable, mrInputVariable);

            KRATOS_INFO_IF("TemporalValueVarianceMethod", this->GetEchoLevel() > 0)
                << "Initialized temporal value variance method for "
                << mrInputVariable.Name() << " input variable with "
                << mrOutputMeanVariable.Name() << " mean variable and "
                << mrOutputVarianceVariable.Name() << " variance variable for "
                << this->GetModelPart().Name() << ".\n";
        }

    private:
        const Variable<TDataType>& mrInputVariable;
        const Variable<TDataType>& mrOutputMeanVariable;
        const Variable<TDataType>& mrOutputVarianceVariable;
    };

    template <class TDataType>
    class NormMethod : public TemporalMethod
    {
    public:
        KRATOS_CLASS_POINTER_DEFINITION(NormMethod);

        NormMethod(
            ModelPart& rModelPart,
            const std::string& rNormType,
            const Variable<TDataType>& rInputVariable,
            const int EchoLevel,
            const Variable<double>& rOutputMeanVariable,
            const Variable<double>& rOutputVarianceVariable)
            : TemporalMethod(rModelPart, EchoLevel),
              mNormType(rNormType),
              mrInputVariable(rInputVariable),
              mrOutputMeanVariable(rOutputMeanVariable),
              mrOutputVarianceVariable(rOutputVarianceVariable)
        {
            KRATOS_TRY

            KRATOS_ERROR_IF(rOutputMeanVariable == rOutputVarianceVariable) << "Same variable is given for mean and variance in norm variance method with input variable "
                                                                            << rInputVariable
                                                                                   .Name()
                                                                            << ". Please provide two different variables. [ variable = "
                                                                            << rOutputMeanVariable
                                                                                   .Name()
                                                                            << " ].\n";

            KRATOS_CATCH("");
        }

        void CalculateStatistics() override
        {
            TContainerType& r_container =
                MethodUtilities::GetDataContainer<TContainerType>(this->GetModelPart());

            const auto& norm_method =
                MethodUtilities::GetNormMethod(mrInputVariable, mNormType);

            const double delta_time = this->GetDeltaTime();
            const double old_total_time = this->GetTotalTime();
            const double total_time = old_total_time + delta_time;

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
                    input_norm_value, delta_time, old_total_time, total_time);
            }

            KRATOS_INFO_IF("TemporalNormVarianceMethod", this->GetEchoLevel() > 1)
                << "Calculated temporal norm variance for " << mrInputVariable.Name()
                << " input variable with " << mrOutputMeanVariable.Name()
                << " mean variable and " << mrOutputVarianceVariable.Name()
                << " variance variable for " << this->GetModelPart().Name() << ".\n";
        }

        // norm output variable initialization
        void InitializeStatisticsVariables() override
        {
            TContainerType& r_container =
                MethodUtilities::GetDataContainer<TContainerType>(this->GetModelPart());

            auto& initializer_method =
                TemporalMethodUtilities::InitializeVariables<TContainerType, TContainerItemType, TDataStorageFunctor>;
            initializer_method(r_container, mrOutputMeanVariable, 0.0);
            initializer_method(r_container, mrOutputVarianceVariable, 0.0);

            KRATOS_INFO_IF("TemporalNormVarianceMethod", this->GetEchoLevel() > 0)
                << "Initialized temporal norm variance method for "
                << mrInputVariable.Name() << " input variable with "
                << mrOutputMeanVariable.Name() << " mean variable and "
                << mrOutputVarianceVariable.Name() << " variance variable for "
                << this->GetModelPart().Name() << ".\n";
        }

    private:
        const std::string mNormType;
        const Variable<TDataType>& mrInputVariable;
        const Variable<double>& mrOutputMeanVariable;
        const Variable<double>& mrOutputVarianceVariable;
    };

    std::vector<TemporalMethod::Pointer> static CreateTemporalMethodObject(
        ModelPart& rModelPart, const std::string& rNormType, const int EchoLevel, Parameters Params)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
            {
                "input_variables"           : [],
                "output_mean_variables"     : [],
                "output_variance_variables" : []
            })");
        Params.RecursivelyValidateAndAssignDefaults(default_parameters);

        const std::vector<std::string>& input_variable_names_list =
            Params["input_variables"].GetStringArray();
        const std::vector<std::string>& output_variable_1_names_list =
            Params["output_mean_variables"].GetStringArray();
        const std::vector<std::string>& output_variable_2_names_list =
            Params["output_variance_variables"].GetStringArray();

        std::vector<TemporalMethod::Pointer> method_list;
        if (rNormType == "none") // for non norm types
        {
            MethodUtilities::CheckInputOutputVariables(
                input_variable_names_list, output_variable_1_names_list);
            MethodUtilities::CheckInputOutputVariables(
                input_variable_names_list, output_variable_2_names_list);

            const int number_of_variables = input_variable_names_list.size();
            for (int i = 0; i < number_of_variables; ++i)
            {
                const std::string& r_variable_input_name = input_variable_names_list[i];
                const std::string& r_variable_1_output_name =
                    output_variable_1_names_list[i];
                const std::string& r_variable_2_output_name =
                    output_variable_2_names_list[i];
                ADD_TEMPORAL_VALUE_METHOD_TWO_OUTPUT_VARIABLE_OBJECT(
                    rModelPart, rNormType, r_variable_input_name, EchoLevel,
                    r_variable_1_output_name, r_variable_2_output_name, method_list, ValueMethod)
            }
        }
        else // for values with norms
        {
            MethodUtilities::CheckVariableType<double>(output_variable_1_names_list);
            MethodUtilities::CheckVariableType<double>(output_variable_2_names_list);

            const int number_of_variables = input_variable_names_list.size();
            for (int i = 0; i < number_of_variables; ++i)
            {
                const std::string& r_variable_input_name = input_variable_names_list[i];
                const std::string& r_variable_1_output_name =
                    output_variable_1_names_list[i];
                const std::string& r_variable_2_output_name =
                    output_variable_2_names_list[i];
                ADD_TEMPORAL_NORM_METHOD_TWO_OUTPUT_VARIABLE_OBJECT(
                    rModelPart, rNormType, r_variable_input_name, EchoLevel,
                    r_variable_1_output_name, r_variable_2_output_name, method_list, NormMethod)
            }
        }

        return method_list;

        KRATOS_CATCH("");
    }

private:
    template <class TDataType>
    void static CalculateMeanAndVariance(
        TDataType& rMean,
        TDataType& rVariance,
        const TDataType& rNewDataPoint,
        const double DeltaTime,
        const double OldTotalTime,
        const double CurrentTotalTime)
    {
        const TDataType new_mean =
            (rMean * OldTotalTime + rNewDataPoint * DeltaTime) * (1.0 / CurrentTotalTime);
        rVariance =
            ((rVariance + MethodUtilities::RaiseToPower<TDataType>(rMean, 2)) * OldTotalTime +
             MethodUtilities::RaiseToPower<TDataType>(rNewDataPoint, 2) * DeltaTime) *
                (1 / CurrentTotalTime) -
            MethodUtilities::RaiseToPower<TDataType>(new_mean, 2);
        rMean = new_mean;
    }
};
} // namespace TemporalMethods
} // namespace Kratos

#endif // KRATOS_TEMPORAL_VARIANCE_METHOD_H_INCLUDED