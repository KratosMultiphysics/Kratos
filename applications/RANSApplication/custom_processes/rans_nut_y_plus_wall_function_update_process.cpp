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

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_nut_y_plus_wall_function_update_process.h"

namespace Kratos
{
RansNutYPlusWallFunctionUpdateProcess::RansNutYPlusWallFunctionUpdateProcess(
    Model& rModel,
    Parameters rParameters)
: mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mEchoLevel = rParameters["echo_level"].GetInt();
    mModelPartName = rParameters["model_part_name"].GetString();
    mMinValue = rParameters["min_value"].GetDouble();

    KRATOS_CATCH("");
}

RansNutYPlusWallFunctionUpdateProcess::RansNutYPlusWallFunctionUpdateProcess(
    Model& rModel,
    const std::string& rModelPartName,
    const double MinValue,
    const int EchoLevel)
    : mrModel(rModel),
      mModelPartName(rModelPartName),
      mMinValue(MinValue),
      mEchoLevel(EchoLevel)
{
}

int RansNutYPlusWallFunctionUpdateProcess::Check()
{
    KRATOS_TRY

    const auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(VISCOSITY))
        << "VISCOSITY is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(TURBULENT_VISCOSITY))
        << "TURBULENT_VISCOSITY is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    return 0;

    KRATOS_CATCH("");
}

void RansNutYPlusWallFunctionUpdateProcess::ExecuteInitialize()
{
    RansCalculationUtilities::CalculateNumberOfNeighbourEntities<ModelPart::ConditionsContainerType>(
        mrModel.GetModelPart(mModelPartName), NUMBER_OF_NEIGHBOUR_CONDITIONS);

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Calculated number of neighbour conditions in " << mModelPartName << ".\n";
}

void RansNutYPlusWallFunctionUpdateProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    if (!mIsInitialized) {
        this->ExecuteAfterCouplingSolveStep();
        mIsInitialized = true;
    }

    KRATOS_CATCH("");
}

void RansNutYPlusWallFunctionUpdateProcess::ExecuteAfterCouplingSolveStep()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);
    VariableUtils().SetHistoricalVariableToZero(TURBULENT_VISCOSITY,
                                                r_model_part.Nodes());

    const double von_karman = r_model_part.GetProcessInfo()[VON_KARMAN];

    block_for_each(r_model_part.Conditions(), [&](ModelPart::ConditionType& rCondition) {
        auto& r_geometry = rCondition.GetGeometry();
        const auto& r_properties = rCondition.GetProperties();
        const double y_plus_limit = r_properties[RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT];
        const double y_plus = std::max(rCondition.GetValue(RANS_Y_PLUS), y_plus_limit);

        auto& r_parent_element = r_geometry.GetValue(NEIGHBOUR_ELEMENTS)[0];
        auto constitutive_law = r_parent_element.GetProperties().GetValue(CONSTITUTIVE_LAW);
        const auto& r_elem_properties = r_parent_element.GetProperties();
        Vector gauss_weights;
        Matrix shape_functions;
        RansCalculationUtilities::CalculateConditionGeometryData(
            r_geometry, GeometryData::IntegrationMethod::GI_GAUSS_1,
            gauss_weights, shape_functions);

        ConstitutiveLaw::Parameters cl_parameters(
            r_geometry, r_elem_properties, r_model_part.GetProcessInfo());
        cl_parameters.SetShapeFunctionsValues(row(shape_functions, 0));

        double nu;
        constitutive_law->CalculateValue(cl_parameters, EFFECTIVE_VISCOSITY, nu);
        nu /= r_elem_properties.GetValue(DENSITY);

        for (IndexType i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node) {
            auto& r_node = r_geometry[i_node];
            r_node.SetLock();
            r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY) +=
                von_karman * y_plus * nu;
            r_node.UnSetLock();
        }
    });

    r_model_part.GetCommunicator().AssembleCurrentData(TURBULENT_VISCOSITY);

    block_for_each(r_model_part.Nodes(), [&](ModelPart::NodeType& rNode) {
        double& r_nut = rNode.FastGetSolutionStepValue(TURBULENT_VISCOSITY);
        const double number_of_neighbour_conditions =
            rNode.GetValue(NUMBER_OF_NEIGHBOUR_CONDITIONS);
        r_nut = RansCalculationUtilities::SoftMax(
            r_nut / number_of_neighbour_conditions, mMinValue);
    });

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 1)
        << "Calculated wall function based nu_t for " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

std::string RansNutYPlusWallFunctionUpdateProcess::Info() const
{
    return std::string("RansNutYPlusWallFunctionUpdateProcess");
}

void RansNutYPlusWallFunctionUpdateProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansNutYPlusWallFunctionUpdateProcess::PrintData(std::ostream& rOStream) const
{
}

const Parameters RansNutYPlusWallFunctionUpdateProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"      : 0,
            "min_value"       : 1e-18
        })");
    return default_parameters;
}

} // namespace Kratos.
