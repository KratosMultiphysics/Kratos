//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <string>

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Include base h
#include "gauss_point_kreisselmeier_aggregation_response_function.h"

namespace Kratos
{

template<class TDataType>
void GaussPointKreisselmeierAggregationResponseFunction::SetGaussPointVariable(
    const Variable<TDataType>*& rpVariable,
    const std::string& rVariableName,
    const std::string& rMsg)
{
    KRATOS_TRY

    if (rVariableName == "") {
        KRATOS_INFO_IF("GaussPointKreisselmeierAggregationResponseFunction", mEchoLevel > 0) << rMsg << std::endl;
    } else {
        rpVariable = &KratosComponents<Variable<TDataType>>::Get(rVariableName);
    }

    KRATOS_CATCH("")
}

template<class TEntityType>
void GaussPointKreisselmeierAggregationResponseFunction::CalculateGaussPointDerivatives(
    Vector& rOutput,
    TEntityType& rEntity,
    const double Factor,
    const Variable<Matrix>* pVariable,
    const Matrix& rResidualGradient,
    const ProcessInfo& rProcessInfo) const
{

    KRATOS_TRY

    if (rOutput.size() != rResidualGradient.size1()) {
        rOutput.resize(rResidualGradient.size1(), false);
    }

    rOutput.clear();

    const auto& p_entity_ks_pre_factor = mKSPrefactors.find(rEntity.Id());
    if (pVariable && p_entity_ks_pre_factor != mKSPrefactors.end()) {
        KRATOS_ERROR_IF_NOT(mAreKSPrefactorsInitialized)
            << "Kreisselmeier aggregation prefactors are not initialized. "
            << "Please run GaussPointKreisselmeierAggregationResponseFunction::CalculateValue to initialize them.\n";

        // may be we can move these matrices to thread local storage
        Matrix gauss_point_derivatives;
        if (mGradientMode == GradientMode::SEMI_ANALITIC) {
            #pragma omp critical
            {
                rEntity.Calculate(*pVariable, gauss_point_derivatives, rProcessInfo);
            }
        } else {
            rEntity.Calculate(*pVariable, gauss_point_derivatives, rProcessInfo);
        }

        for (IndexType i = 0; i < gauss_point_derivatives.size2(); ++i) {
            noalias(rOutput) += column(gauss_point_derivatives, i);
        }

        noalias(rOutput) = rOutput * (Factor * p_entity_ks_pre_factor->second / (mSumKSPrefactors * gauss_point_derivatives.size2() * mGaussPointValueScalingFactor));
    }

    KRATOS_CATCH("");

}

//template instantiations
template void GaussPointKreisselmeierAggregationResponseFunction::SetGaussPointVariable<double>(const Variable<double>*&, const std::string&, const std::string&);
template void GaussPointKreisselmeierAggregationResponseFunction::SetGaussPointVariable<Vector>(const Variable<Vector>*&, const std::string&, const std::string&);

template void GaussPointKreisselmeierAggregationResponseFunction::CalculateGaussPointDerivatives<ModelPart::ElementType>(Vector&, ModelPart::ElementType&, const double, const Variable<Matrix>*,  const Matrix&, const ProcessInfo&) const;

GaussPointKreisselmeierAggregationResponseFunction::GaussPointKreisselmeierAggregationResponseFunction(
    Parameters Settings,
    ModelPart& rModelPart)
    : mrModelPart(rModelPart)
{
    KRATOS_TRY;

    Parameters default_settings(R"(
    {
        "gauss_point_values_scalar_variable"            : "PLEASE_SPECIFY_A_SCALAR_VARIABLE",
        "gauss_point_gradient_matrix_variable"          : "PLEASE_SPECIFY_A_MATRIX_VARIABLE",
        "gauss_point_first_derivatives_matrix_variable" : "PLEASE_SPECIFY_A_MATRIX_VARIABLE",
        "gauss_point_second_derivatives_matrix_variable": "PLEASE_SPECIFY_A_MATRIX_VARIABLE",
        "gauss_point_shape_derivatives_matrix_variable" : "PLEASE_SPECIFY_A_MATRIX_VARIABLE",
        "gauss_point_value_scaling_factor"              : "max",
        "aggregation_penalty"                           : 50.0,
        "echo_level"                                    : 0,
        "gradient_mode"                                 : "semi_analytic",
        "gradient_mode_settings": {
            "perturbation_variable_name"           : "PLEASE_SPECIFY_PERTURBATION_SCALAR_VARIABLE_NAME",
            "perturbation_size"                    : 1e-6,
            "design_variable_name_storage_variable": "PLEASE_SPECIFY_STRING_VARIABLE"
        }
    })");

    if (Settings.Has("gauss_point_value_scaling_factor") && Settings["gauss_point_value_scaling_factor"].IsDouble()) {
        default_settings["gauss_point_value_scaling_factor"].SetDouble(1.0);
    }
    Settings.ValidateAndAssignDefaults(default_settings);

    mAggregationPenalty = Settings["aggregation_penalty"].GetDouble();

    // set the gradient mode settings
    const std::string& gradient_mode = Settings["gradient_mode"].GetString();
    if (gradient_mode == "semi_analytic") {
        Parameters default_semi_analytic_gradient_mode_settings(R"(
        {
            "perturbation_variable_name"           : "PLEASE_SPECIFY_PERTURBATION_SCALAR_VARIABLE_NAME",
            "perturbation_size"                    : 1e-6,
            "design_variable_name_storage_variable": "PLEASE_SPECIFY_STRING_VARIABLE"
        })");

        Settings["gradient_mode_settings"].ValidateAndAssignDefaults(default_semi_analytic_gradient_mode_settings);

        mGradientMode = GradientMode::SEMI_ANALITIC;
        mpPerturbationVariable = &KratosComponents<Variable<double>>::Get(Settings["gradient_mode_settings"]["perturbation_variable_name"].GetString());
        mpDeisgnVariableNameStorageVariable = &KratosComponents<Variable<std::string>>::Get(Settings["gradient_mode_settings"]["design_variable_name_storage_variable"].GetString());
        mPerturbationSize = Settings["gradient_mode_settings"]["perturbation_size"].GetDouble();
    } else if (gradient_mode == "analytic") {
        mGradientMode = GradientMode::ANALYTIC;
    } else {
        KRATOS_ERROR << "Unsupported gradient mode requested. [ gradient_mode " << gradient_mode << ". Supported gradient modes:"
                     << "\n\tsemi_analytic"
                     << "\n\tanalytic" << std::endl;
    }

    // set variables
    mpGaussPointValueScalarVariable = &KratosComponents<Variable<double>>::Get(Settings["gauss_point_values_scalar_variable"].GetString());
    mpGaussPointValueShapeDerivativeVariable = &KratosComponents<Variable<Matrix>>::Get(Settings["gauss_point_shape_derivatives_matrix_variable"].GetString());

    // rest of the variables can be optional as well depending on the type of problem being solved.
    SetGaussPointVariable(mpGaussPointValueGradientVariable, Settings["gauss_point_gradient_matrix_variable"].GetString(), "No gauss_point_gradient_matrix_variable provided.");
    SetGaussPointVariable(mpGaussPointValueFirstDerivativeVariable, Settings["gauss_point_first_derivatives_matrix_variable"].GetString(), "No gauss_point_first_derivatives_matrix_variable provided.");
    SetGaussPointVariable(mpGaussPointValueSecondDerivativeVariable, Settings["gauss_point_second_derivatives_matrix_variable"].GetString(), "No gauss_point_second_derivatives_matrix_variable provided.");

    // set gauss point scaling factor
    if (Settings["gauss_point_value_scaling_factor"].IsDouble()) {
        mGaussPointValueScalingFactor = Settings["gauss_point_value_scaling_factor"].GetDouble();
        KRATOS_ERROR_IF(mGaussPointValueScalingFactor <= 0.0)
            << "\"gauss_point_value_scaling_factor\" needs to be a positive value. [ gauss_point_value_scaling_factor = "
            << mGaussPointValueScalingFactor << " ].\n";
    }

    KRATOS_CATCH("");
}

void GaussPointKreisselmeierAggregationResponseFunction::Initialize()
{
    KRATOS_TRY

    if (mGradientMode == GradientMode::SEMI_ANALITIC) {
        mrModelPart.GetProcessInfo()[*mpPerturbationVariable] = mPerturbationSize;
        KRATOS_INFO_IF("GaussPointKreisselmeierAggregationResponseFunction", mEchoLevel > 0)
                << "Initialized " << mpPerturbationVariable->Name() << " with " << mPerturbationSize
                << " in " << mrModelPart.FullName()
                << " process info for semi-analytic derivative computation.\n";
    }

    KRATOS_CATCH("");
}

double GaussPointKreisselmeierAggregationResponseFunction::CalculateValue(ModelPart& rModelPart)
{
    KRATOS_TRY

    const auto& r_process_info = mrModelPart.GetProcessInfo();

    for (const auto& r_element : mrModelPart.Elements()) {
        mKSPrefactors[r_element.Id()] = 0.0;
    }

    const double max_mean_gp_value = block_for_each<MaxReduction<double>>(mrModelPart.Elements(), std::vector<double>(), [&](auto& rElement, std::vector<double>& rTLS) -> double {
        rElement.CalculateOnIntegrationPoints(*mpGaussPointValueScalarVariable, rTLS, r_process_info);

        double mean_gp_value = 0.0;
        for (IndexType i = 0; i < rTLS.size(); ++i) {
            mean_gp_value += rTLS[i];
        }

        mean_gp_value /= rTLS.size();
        mKSPrefactors.find(rElement.Id())->second = mean_gp_value;

        return mean_gp_value;
    });

    if (mGaussPointValueScalingFactor == 0.0) {
        mGaussPointValueScalingFactor = max_mean_gp_value;
        KRATOS_INFO_IF("GaussPointKreisselmeierAggregationResponseFunction", mEchoLevel > 0)
            << "Using maximum gauss point value as gauss_point_scaling_factor [ gauss_point_scaling_factor = "
            << mGaussPointValueScalingFactor << " ].\n";
    }

    mSumKSPrefactors = block_for_each<SumReduction<double>>(mrModelPart.Elements(), [&](const auto& rElement) -> double {
        double& r_element_ks_prefactor = mKSPrefactors.find(rElement.Id())->second;
        r_element_ks_prefactor = std::exp(mAggregationPenalty * r_element_ks_prefactor / mGaussPointValueScalingFactor);
        return r_element_ks_prefactor;
    });

    mAreKSPrefactorsInitialized = true;

    return std::log(mSumKSPrefactors) / mAggregationPenalty;

    KRATOS_CATCH("");
}

void GaussPointKreisselmeierAggregationResponseFunction::CalculateGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    CalculateGaussPointDerivatives(rResponseGradient, const_cast<ModelPart::ElementType&>(rAdjointElement), -1.0 ,mpGaussPointValueGradientVariable, rResidualGradient, rProcessInfo);

    KRATOS_CATCH("");
}



void GaussPointKreisselmeierAggregationResponseFunction::CalculateGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    if(rResponseGradient.size() != rResidualGradient.size1()) {
        rResponseGradient.resize(rResidualGradient.size1(), false);
    }
    rResponseGradient.clear();

    KRATOS_CATCH("");
}

void GaussPointKreisselmeierAggregationResponseFunction::CalculateFirstDerivativesGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    CalculateGaussPointDerivatives(rResponseGradient, const_cast<ModelPart::ElementType&>(rAdjointElement), -1.0, mpGaussPointValueFirstDerivativeVariable, rResidualGradient, rProcessInfo);

    KRATOS_CATCH("");
}

void GaussPointKreisselmeierAggregationResponseFunction::CalculateFirstDerivativesGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    if(rResponseGradient.size() != rResidualGradient.size1()) {
        rResponseGradient.resize(rResidualGradient.size1(), false);
    }

    rResponseGradient.clear();

    KRATOS_CATCH("");
}

void GaussPointKreisselmeierAggregationResponseFunction::CalculateSecondDerivativesGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    CalculateGaussPointDerivatives(rResponseGradient, const_cast<ModelPart::ElementType&>(rAdjointElement), -1.0, mpGaussPointValueSecondDerivativeVariable, rResidualGradient, rProcessInfo);

    KRATOS_CATCH("");
}

void GaussPointKreisselmeierAggregationResponseFunction::CalculateSecondDerivativesGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    if(rResponseGradient.size() != rResidualGradient.size1()) {
        rResponseGradient.resize(rResidualGradient.size1(), false);
    }

    rResponseGradient.clear();

    KRATOS_CATCH("");
}

void GaussPointKreisselmeierAggregationResponseFunction::CalculatePartialSensitivity(
    Element& rAdjointElement,
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    if (mGradientMode == GradientMode::SEMI_ANALITIC) {
        rAdjointElement.SetValue(*mpDeisgnVariableNameStorageVariable, rVariable.Name());
    }

    if (rVariable == SHAPE_SENSITIVITY) {
        CalculateGaussPointDerivatives(rSensitivityGradient, rAdjointElement, 1.0, mpGaussPointValueShapeDerivativeVariable, rSensitivityMatrix, rProcessInfo);
        KRATOS_ERROR_IF(rSensitivityGradient.size() != rSensitivityMatrix.size1()) << "Size of partial stress design variable derivative does not fit!" << std::endl;
    } else {
        KRATOS_ERROR << "CalculatePartialSensitivity for " << rVariable.Name() << " is not defined.\n";
    }

    if (mGradientMode == GradientMode::SEMI_ANALITIC) {
        rAdjointElement.SetValue(*mpDeisgnVariableNameStorageVariable, "");
    }

    KRATOS_CATCH("");
}

void GaussPointKreisselmeierAggregationResponseFunction::CalculatePartialSensitivity(
    Condition& rAdjointCondition,
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    if (rVariable == SHAPE_SENSITIVITY) {
        if(rSensitivityGradient.size() != rSensitivityMatrix.size1()) {
            rSensitivityGradient.resize(rSensitivityMatrix.size1(), false);
        }
        rSensitivityGradient.clear();
    } else {
        KRATOS_ERROR << "CalculatePartialSensitivity for " << rVariable.Name() << " is not defined.\n";
    }

    KRATOS_CATCH("");
}

} /* namespace Kratos.*/

