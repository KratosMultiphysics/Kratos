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

GaussPointKreisselmeierAggregationResponseFunction::GaussPointKreisselmeierAggregationResponseFunction(
    Parameters Settings,
    ModelPart& rModelPart)
    : mrModelPart(rModelPart)
{
    KRATOS_TRY;

    Parameters default_settings(R"(
    {
        "model_part_name"                   : "PLEASE_SPECIFY_VOLUME_MODEL_PART_NAME",
        "gauss_point_values_scalar_variable": "PLEASE_SPECIFY_A_SCALAR_VARIABLE",
        "rho"                               : 1.0,
        "gauss_point_value_scaling_factor"  : 0.0,
        "echo_level"                        : 0,
        "gradient_mode"                     : "semi_analytic",
        "gradient_mode_settings"            : {
            "perturbation_variable_name"                             : "PLEASE_SPECIFY_PERTURBATION_SCALAR_VARIABLE_NAME",
            "perturbation_size"                                      : 1e-6,
            "list_of_scalar_primal_state_variables"                  : ["PLEASE_SPECIFY_PRIMAL_STATE_SCALAR_VARIABLES_LIST"],
            "list_of_scalar_primal_state_first_derivative_variables" : [],
            "list_of_scalar_primal_state_second_derivative_variables": []
        }
    })");

    Settings.RecursivelyValidateAndAssignDefaults(default_settings);

    mCriticalModelPartName = Settings["model_part_name"].GetString();
    mPerturbationVariableName = Settings["gradient_mode_settings"]["perturbation_variable_name"].GetString();
    KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(mPerturbationVariableName))
        << mPerturbationVariableName << " is not found in the scalar variables list for perturbation_variable_name.\n";

    mGaussPointValueScalingFactor = Settings["gauss_point_value_scaling_factor"].GetDouble();
    mRho = Settings["rho"].GetDouble();
    mStepSize = Settings["gradient_mode_settings"]["perturbation_size"].GetDouble();

    const auto& r_gauss_point_values_scalar_variable = Settings["gauss_point_values_scalar_variable"].GetString();
    KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(r_gauss_point_values_scalar_variable))
        << r_gauss_point_values_scalar_variable << " is not found in the scalar variables list for gauss_point_values_scalar_variable.\n";
    mpGaussPointValueScalarVariable = &(KratosComponents<Variable<double>>::Get(r_gauss_point_values_scalar_variable));

    for (const auto& r_scalar_primal_state_variable_name : Settings["gradient_mode_settings"]["list_of_scalar_primal_state_variables"].GetStringArray()) {
        KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(r_scalar_primal_state_variable_name))
            << r_scalar_primal_state_variable_name << " is not found in the scalar variables list used in list_of_scalar_primal_state_variables.\n";
        mPrimalStateScalarVariablePointersList.push_back(&(KratosComponents<Variable<double>>::Get(r_scalar_primal_state_variable_name)));
    }

    for (const auto& r_scalar_primal_state_variable_name : Settings["gradient_mode_settings"]["list_of_scalar_primal_state_first_derivative_variables"].GetStringArray()) {
        KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(r_scalar_primal_state_variable_name))
            << r_scalar_primal_state_variable_name << " is not found in the scalar variables list used in list_of_scalar_primal_state_first_derivative_variables.\n";
        mPrimalStateFirstDerivativeScalarVariablePointersList.push_back(&(KratosComponents<Variable<double>>::Get(r_scalar_primal_state_variable_name)));
    }

    for (const auto& r_scalar_primal_state_variable_name : Settings["gradient_mode_settings"]["list_of_scalar_primal_state_second_derivative_variables"].GetStringArray()) {
        KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(r_scalar_primal_state_variable_name))
            << r_scalar_primal_state_variable_name << " is not found in the scalar variables list used in list_of_scalar_primal_state_second_derivative_variables.\n";
        mPrimalStateSecondDerivativeScalarVariablePointersList.push_back(&(KratosComponents<Variable<double>>::Get(r_scalar_primal_state_variable_name)));
    }

    KRATOS_CATCH("");
}

void GaussPointKreisselmeierAggregationResponseFunction::Initialize()
{
    KRATOS_TRY

    const auto& r_perturbation_variable = KratosComponents<Variable<double>>::Get(mPerturbationVariableName);
    mrModelPart.GetProcessInfo()[r_perturbation_variable] = mStepSize;

    KRATOS_CATCH("");
}

double GaussPointKreisselmeierAggregationResponseFunction::CalculateValue(ModelPart& rModelPart)
{
    KRATOS_TRY

    const auto& r_process_info = rModelPart.GetProcessInfo();
    auto& r_critical_model_part = rModelPart.GetSubModelPart(mCriticalModelPartName);

    for (const auto& r_element : r_critical_model_part.Elements()) {
        mKSPrefactors[r_element.Id()] = 0.0;
    }

    const double max_mean_gp_value = block_for_each<MaxReduction<double>>(r_critical_model_part.Elements(), std::vector<double>(), [&](ModelPart::ElementType& rElement, std::vector<double>& rTLS) -> double {
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
    }

    mSumKSPrefactors = block_for_each<SumReduction<double>>(r_critical_model_part.Elements(), [&](ModelPart::ElementType& rElement) -> double {
        double& r_element_ks_prefactor = mKSPrefactors.find(rElement.Id())->second;
        r_element_ks_prefactor = std::exp(mRho * r_element_ks_prefactor / mGaussPointValueScalingFactor);
        return r_element_ks_prefactor;
    });

    mAreKSPrefactorsInitialized = true;

    return std::log(mSumKSPrefactors) / mRho;

    KRATOS_CATCH("");
}

void GaussPointKreisselmeierAggregationResponseFunction::CalculateFiniteDifferenceStateVariableSensitivities(
    Vector& rOutput,
    ModelPart::ElementType& rElement,
    const std::vector<const Variable<double>*>& rDerivativeVariablePointersList,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY
    auto& r_geometry = rElement.GetGeometry();

    const IndexType number_of_nodes = r_geometry.PointsNumber();
    const IndexType number_of_dofs_per_node = rDerivativeVariablePointersList.size();
    const IndexType local_derivative_size = number_of_nodes * number_of_dofs_per_node;

    std::vector<double> ref_gp_values;
    rElement.CalculateOnIntegrationPoints(*mpGaussPointValueScalarVariable, ref_gp_values, rProcessInfo);

    if (rOutput.size() != local_derivative_size) {
        rOutput.resize(local_derivative_size, false);
    }

    rOutput.clear();

    #pragma omp critical
    {
        std::vector<double> gp_values;
        for (IndexType i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node) {

            auto& r_node = r_geometry[i_node];
            for (IndexType i_var = 0; i_var < number_of_dofs_per_node; ++i_var) {
                const auto p_variable = rDerivativeVariablePointersList[i_var];
                double& variable_value = r_node.FastGetSolutionStepValue(*p_variable);

                variable_value += mStepSize;

                rElement.CalculateOnIntegrationPoints(*mpGaussPointValueScalarVariable, gp_values, rProcessInfo);
                for (IndexType i = 0; i < gp_values.size(); ++i) {
                    rOutput(i_node * number_of_dofs_per_node + i_var) += (gp_values[i] - ref_gp_values[i]) / mStepSize;
                }

                rOutput(i_node * number_of_dofs_per_node + i_var) /= ref_gp_values.size();

                variable_value -= mStepSize;
            }
        }
    }

    KRATOS_CATCH("");
}

void GaussPointKreisselmeierAggregationResponseFunction::CalculateFiniteDifferenceShapeVariableSensitivities(
    Vector& rOutput,
    ModelPart::ElementType& rElement,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY
    auto& r_geometry = rElement.GetGeometry();

    const IndexType number_of_nodes = r_geometry.PointsNumber();
    const IndexType number_of_dofs_per_node = rProcessInfo[DOMAIN_SIZE];
    const IndexType local_derivative_size = number_of_nodes * number_of_dofs_per_node;

    std::vector<double> ref_gp_values;
    rElement.CalculateOnIntegrationPoints(*mpGaussPointValueScalarVariable, ref_gp_values, rProcessInfo);

    if (rOutput.size() != local_derivative_size) {
        rOutput.resize(local_derivative_size, false);
    }

    rOutput.clear();
    #pragma omp critical
    {
        std::vector<double> gp_values;
        for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {

            auto& r_node = r_geometry[i_node];
            for (IndexType i_var = 0; i_var < number_of_dofs_per_node; ++i_var) {
                double& position = r_node.Coordinates()[i_var];
                double& initial_position = r_node.GetInitialPosition()[i_var];

                position += mStepSize;
                initial_position += mStepSize;

                rElement.CalculateOnIntegrationPoints(*mpGaussPointValueScalarVariable, gp_values, rProcessInfo);
                for (IndexType i = 0; i < gp_values.size(); ++i) {
                    rOutput(i_node * number_of_dofs_per_node + i_var) += (gp_values[i] - ref_gp_values[i]) / mStepSize;
                }

                position -= mStepSize;
                initial_position -= mStepSize;
            }
        }
    }

    KRATOS_CATCH("");
}

void GaussPointKreisselmeierAggregationResponseFunction::CalculateStateDerivative(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const std::vector<const Variable<double>*>& rDerivativeVariablesList,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mAreKSPrefactorsInitialized)
        << "GaussPointKreisselmeierAggregationResponseFunction::CalculateGradient: Prefactors missing. First calculate value before calculating gradients!"
        << std::endl;

    if (mKSPrefactors.find(rAdjointElement.Id()) != mKSPrefactors.end() && rDerivativeVariablesList.size() != 0) {
        ModelPart::ElementType& r_element = *(mrModelPart.pGetElement(rAdjointElement.Id()));
        CalculateFiniteDifferenceStateVariableSensitivities(rResponseGradient, r_element, rDerivativeVariablesList, rProcessInfo);
        rResponseGradient *= mKSPrefactors.find(r_element.Id())->second / (-mGaussPointValueScalingFactor * mSumKSPrefactors);
    } else {
        if(rResponseGradient.size() != rResidualGradient.size1()) {
            rResponseGradient.resize(rResidualGradient.size1(), false);
        }
        rResponseGradient.clear();
    }

    KRATOS_CATCH("");
}

void GaussPointKreisselmeierAggregationResponseFunction::CalculateGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    CalculateStateDerivative(rAdjointElement, rResidualGradient, rResponseGradient, mPrimalStateScalarVariablePointersList, rProcessInfo);

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

    CalculateStateDerivative(rAdjointElement, rResidualGradient, rResponseGradient, mPrimalStateFirstDerivativeScalarVariablePointersList, rProcessInfo);

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

    CalculateStateDerivative(rAdjointElement, rResidualGradient, rResponseGradient, mPrimalStateSecondDerivativeScalarVariablePointersList, rProcessInfo);

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

    if (rVariable == SHAPE_SENSITIVITY) {
        if (mKSPrefactors.find(rAdjointElement.Id()) != mKSPrefactors.end()) {
            CalculateFiniteDifferenceShapeVariableSensitivities(rSensitivityGradient, rAdjointElement, rProcessInfo);
            rSensitivityGradient *= mKSPrefactors.find(rAdjointElement.Id())->second / (-mGaussPointValueScalingFactor * mSumKSPrefactors);
        } else {
            if(rSensitivityGradient.size() != rSensitivityMatrix.size1()) {
                rSensitivityGradient.resize(rSensitivityMatrix.size1(), false);
            }
            rSensitivityGradient.clear();
        }
    } else {
        KRATOS_ERROR << "CalculatePartialSensitivity for " << rVariable.Name() << " is not defined.\n";
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

