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
#include <functional>

// External includes

// Project includes
#include "containers/model.h"
#include "containers/variable.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "utilities/parallel_utilities.h"

// Application incldues

// Include base h
#include "calculate_sensitivity_metric_process.h"

namespace Kratos
{
CalculateSensitivityMetricProcess::CalculateSensitivityMetricProcess(
    Model& rModel,
    Parameters rParameters)
    : mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mModelPartName = rParameters["model_part_name"].GetString();
    mSensitivityVariableName = rParameters["sensitivity_variable_name"].GetString();
    mIsHistoricalVariable = rParameters["is_historical_variable"].GetBool();

    KRATOS_CATCH("");
}

void CalculateSensitivityMetricProcess::Execute()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);
    const auto& r_process_info = r_model_part.GetProcessInfo();
    const IndexType domain_size = r_process_info[DOMAIN_SIZE];

    const auto& r_sensitivity_variable = KratosComponents<Variable<array_1d<double, 3>>>::Get(mSensitivityVariableName);

    std::function<Array3D(const NodeType& rNode)> sensitivity_getter_method;
    if (mIsHistoricalVariable) {
        KRATOS_ERROR_IF(r_model_part.HasNodalSolutionStepVariable(r_sensitivity_variable))
            << mSensitivityVariableName << " not found in solution step variables list of "
            << mModelPartName << ".\n";

        sensitivity_getter_method = [&](const NodeType& rNode) -> Array3D {
            return rNode.FastGetSolutionStepValue(r_sensitivity_variable);
        };
    } else {
        sensitivity_getter_method = [&](const NodeType& rNode) -> Array3D {
            return rNode.GetValue(r_sensitivity_variable);
        };
    }

    std::function<void(NodeType&, const Array3D&)> metric_setter_method;
    if (domain_size == 2) {
        const auto& r_metric_variable = KratosComponents<Variable<array_1d<double, 3>>>::Get("METRIC_TENSOR_2D");
        metric_setter_method = [&](NodeType& rNode, const Array3D& rValue) {
            rNode.SetValue(r_metric_variable,
                           array_1d<double, 3>({rValue[0], rValue[1], 0.0}));
        };
    } else if (domain_size == 3) {
        const auto& r_metric_variable = KratosComponents<Variable<array_1d<double, 6>>>::Get("METRIC_TENSOR_3D");
        metric_setter_method = [&](NodeType& rNode, const Array3D& rValue) {
            rNode.SetValue(
                r_metric_variable,
                array_1d<double, 6>({rValue[0], rValue[1], rValue[2], 0.0, 0.0, 0.0}));
        };
    } else {
        KRATOS_ERROR << "Unsupported domain size provided. [ domain_size = " << domain_size
                     << " ].\n";
    }

    block_for_each(r_model_part.Nodes(), [&](NodeType& rNode) {
        metric_setter_method(rNode, sensitivity_getter_method(rNode));
    });

    KRATOS_CATCH("");
}

const Parameters CalculateSensitivityMetricProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "model_part_name"          : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "sensitivity_variable_name": "SHAPE_SENSITIVITY",
            "is_historical_variable"   : true
        })");

    return default_parameters;
}

std::string CalculateSensitivityMetricProcess::Info() const
{
    return std::string("CalculateSensitivityMetricProcess");
}

/// Print information about this object.
void CalculateSensitivityMetricProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

/// Print object's data.
void CalculateSensitivityMetricProcess::PrintData(std::ostream& rOStream) const
{

}

} // namespace Kratos.

