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

#if !defined(KRATOS_TEMPORAL_MEAN_METHOD_H_INCLUDED)
#define KRATOS_TEMPORAL_MEAN_METHOD_H_INCLUDED

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
class TemporalMeanMethod
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
            const Variable<TDataType>& rOutputVariable)
            : TemporalMethod(rModelPart, EchoLevel),
              mrInputVariable(rInputVariable),
              mrOutputVariable(rOutputVariable)
        {
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
                TDataType& r_output_value =
                    TDataStorageFunctor<TContainerItemType>()(r_item, mrOutputVariable);
                MethodUtilities::DataTypeSizeChecker(r_input_value, r_output_value);

                TemporalMeanMethod::CalculateMean<TDataType>(
                    r_output_value, r_input_value, delta_time, old_total_time, total_time);
            }

            KRATOS_INFO_IF("TemporalValueMeanMethod", this->GetEchoLevel() > 1)
                << "Calculated temporal value mean for " << mrInputVariable.Name()
                << " input variable with " << mrOutputVariable.Name()
                << " mean variable for " << this->GetModelPart().Name() << ".\n";
        }

        void InitializeStatisticsVariables() override
        {
            TContainerType& r_container =
                MethodUtilities::GetDataContainer<TContainerType>(this->GetModelPart());

            auto& initializer_method =
                TemporalMethodUtilities::InitializeVariables<TContainerType, TContainerItemType, TDataRetrievalFunctor, TDataStorageFunctor, TDataType>;
            initializer_method(r_container, mrOutputVariable, mrInputVariable);

            KRATOS_INFO_IF("TemporalValueMeanMethod", this->GetEchoLevel() > 0)
                << "Initialized temporal value mean method for "
                << mrInputVariable.Name() << " input variable with "
                << mrOutputVariable.Name() << " mean variable for "
                << this->GetModelPart().Name() << ".\n";
        }

    private:
        const Variable<TDataType>& mrInputVariable;
        const Variable<TDataType>& mrOutputVariable;
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
            const Variable<double>& rOutputVariable)
            : TemporalMethod(rModelPart, EchoLevel),
              mNormType(rNormType),
              mrInputVariable(rInputVariable),
              mrOutputVariable(rOutputVariable)
        {
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
                double& r_output_value =
                    TDataStorageFunctor<TContainerItemType>()(r_item, mrOutputVariable);

                TemporalMeanMethod::CalculateMean<double>(
                    r_output_value, input_norm_value, delta_time, old_total_time, total_time);
            }

            KRATOS_INFO_IF("TemporalNormMeanMethod", this->GetEchoLevel() > 1)
                << "Calculated temporal norm mean for " << mrInputVariable.Name()
                << " input variable with " << mrOutputVariable.Name()
                << " mean variable for " << this->GetModelPart().Name() << ".\n";
        }

        void InitializeStatisticsVariables() override
        {
            TContainerType& r_container =
                MethodUtilities::GetDataContainer<TContainerType>(this->GetModelPart());

            auto& initializer_method =
                TemporalMethodUtilities::InitializeVariables<TContainerType, TContainerItemType, TDataStorageFunctor>;
            initializer_method(r_container, mrOutputVariable, 0.0);

            KRATOS_INFO_IF("TemporalNormMeanMethod", this->GetEchoLevel() > 0)
                << "Initialized temporal norm mean method for "
                << mrInputVariable.Name() << " input variable with "
                << mrOutputVariable.Name() << " mean variable for "
                << this->GetModelPart().Name() << ".\n";
        }

    private:
        const std::string mNormType;
        const Variable<TDataType>& mrInputVariable;
        const Variable<double>& mrOutputVariable;
    };

    std::vector<TemporalMethod::Pointer> static CreateTemporalMethodObject(
        ModelPart& rModelPart, const std::string& rNormType, const int EchoLevel, Parameters Params)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
            {
                "input_variables"  : [],
                "output_variables" : []
            })");
        Params.RecursivelyValidateAndAssignDefaults(default_parameters);

        const std::vector<std::string>& input_variable_names_list =
            Params["input_variables"].GetStringArray();
        const std::vector<std::string>& output_variable_names_list =
            Params["output_variables"].GetStringArray();

        std::vector<TemporalMethod::Pointer> method_list;
        if (rNormType == "none") // for non norm types
        {
            MethodUtilities::CheckInputOutputVariables(
                input_variable_names_list, output_variable_names_list);
            const int number_of_variables = input_variable_names_list.size();
            for (int i = 0; i < number_of_variables; ++i)
            {
                const std::string& r_variable_input_name = input_variable_names_list[i];
                const std::string& r_variable_output_name =
                    output_variable_names_list[i];
                ADD_TEMPORAL_VALUE_METHOD_ONE_OUTPUT_VARIABLE_OBJECT(
                    rModelPart, rNormType, r_variable_input_name, EchoLevel,
                    r_variable_output_name, method_list, ValueMethod)
            }
        }
        else // for values with norms
        {
            MethodUtilities::CheckVariableType<double>(output_variable_names_list);

            const int number_of_variables = input_variable_names_list.size();
            for (int i = 0; i < number_of_variables; ++i)
            {
                const std::string& r_variable_input_name = input_variable_names_list[i];
                const std::string& r_variable_output_name =
                    output_variable_names_list[i];
                ADD_TEMPORAL_NORM_METHOD_ONE_OUTPUT_VARIABLE_OBJECT(
                    rModelPart, rNormType, r_variable_input_name, EchoLevel,
                    r_variable_output_name, method_list, NormMethod)
            }
        }

        return method_list;

        KRATOS_CATCH("");
    }

private:
    template <class TDataType>
    void static CalculateMean(
        TDataType& rMean,
        const TDataType& rNewDataPoint,
        const double DeltaTime,
        const double OldTotalTime,
        const double CurrentTotalTime)
    {
        rMean = (rMean * OldTotalTime + rNewDataPoint * DeltaTime) * (1.0 / CurrentTotalTime);
    }
};
} // namespace TemporalMethods
} // namespace Kratos

#endif // KRATOS_TEMPORAL_MEAN_METHOD_H_INCLUDED