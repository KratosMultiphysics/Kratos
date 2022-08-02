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
#include <cmath>

// External includes

// Project includes
#include "includes/checks.h"
#include "utilities/normal_calculation_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"
#include "includes/variables.h"

// Application includes
#include "fluid_dynamics_application_variables.h"

// Include base h
#include "compute_y_plus_process.h"

namespace Kratos
{
const Parameters ComputeYPlusProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "model_part_name"                  : "PLEASE_SPECIFY_MAIN_MODEL_PART_NAME",
            "output_variable_name"             : "Y_PLUS",
            "output_to_elements"               : false,
            "calculate_normals_every_time_step": false,
            "echo_level"                       : 0
        })");

    return default_parameters;
}

ComputeYPlusProcess::ComputeYPlusProcess(
    Model& rModel,
    Parameters rParameters)
    : mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mModelPartName = rParameters["model_part_name"].GetString();
    mOutputVariableName = rParameters["output_variable_name"].GetString();
    mIsOutputStoredInElements = rParameters["output_to_elements"].GetBool();
    mIsCalculatedEveryTimeStep = rParameters["calculate_normals_every_time_step"].GetBool();
    mEchoLevel = rParameters["echo_level"].GetInt();

    mIsNormalsCalculated = false;

    KRATOS_CATCH("");
}

void ComputeYPlusProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    if (!mIsNormalsCalculated || mIsCalculatedEveryTimeStep) {
        NormalCalculationUtils().CalculateNormals<ConditionsContainerType>(r_model_part);

        KRATOS_INFO_IF(this->Info(), !mIsNormalsCalculated || mEchoLevel > 0)
                            << "Calculated normals in " << r_model_part.FullName() << ".\n";
    }

    mIsNormalsCalculated = true;

    KRATOS_CATCH("");
}

void ComputeYPlusProcess::ExecuteFinalizeSolutionStep()
{
    KRATOS_TRY

    CalculateConditionReactions();

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    auto& r_output_variable = KratosComponents<Variable<double>>::Get(mOutputVariableName);

    block_for_each(r_model_part.Conditions(), [&](ConditionType& rCondition) {
        // calculate condition unit normal
        array_1d<double, 3> r_condition_unit_normal = rCondition.GetValue(NORMAL);
        r_condition_unit_normal /= norm_2(r_condition_unit_normal);

        // get parent element for which this condition belongs to
        const auto& r_parent_element = rCondition.GetValue(NEIGHBOUR_ELEMENTS)[0];

        // calculate wall height
        const auto& r_parent_geometry = r_parent_element.GetGeometry();
        const auto& r_condition_geometry = rCondition.GetGeometry();

        const auto& parent_center = r_parent_geometry.Center();
        const auto& condition_center = r_condition_geometry.Center();

        const double y = inner_prod(condition_center - parent_center, r_condition_unit_normal);

        // get fluid properties from parent element
        const auto& r_properties = r_parent_element.GetProperties();
        const double density = r_properties[DENSITY];
        const double kinmeatic_viscosity = r_properties[DYNAMIC_VISCOSITY] / density;

        // calculate reaction
        const array_1d<double, 3>& r_reaction = rCondition.GetValue(REACTION);
        const array_1d<double, 3>& perpendicular_reaction = r_condition_unit_normal * inner_prod(r_reaction, r_condition_unit_normal);
        const array_1d<double, 3>& tangential_reaction = r_reaction - perpendicular_reaction;

        const double shear_stress = norm_2(tangential_reaction) / r_condition_geometry.Area();
        const double u_tau = std::sqrt(shear_stress / density);

        rCondition.SetValue(r_output_variable, u_tau * y / kinmeatic_viscosity);
    });

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)<< "Calculated y_plus in " << r_model_part.FullName() << ".\n";

    KRATOS_CATCH("");
}

int ComputeYPlusProcess::Check()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    KRATOS_ERROR_IF_NOT(r_model_part.HasNodalSolutionStepVariable(REACTION))
                    << "REACTION variable is not found in nodal solution step variables list of "
                    << r_model_part.FullName() << ".\n";

    KRATOS_ERROR_IF_NOT(r_model_part.HasNodalSolutionStepVariable(NORMAL))
                    << "NORMAL variable is not found in nodal solution step variables list of "
                    << r_model_part.FullName() << ".\n";

    block_for_each(r_model_part.Conditions(), [](const ConditionType& rCondition) {
        KRATOS_ERROR_IF_NOT(rCondition.Has(NEIGHBOUR_ELEMENTS)) << "NEIGHBOUR_ELEMENTS is not present in condition with id " << rCondition.Id() << ".\n";
        KRATOS_ERROR_IF_NOT(rCondition.GetValue(NEIGHBOUR_ELEMENTS).size() == 1) << "NEIGHBOUR_ELEMENTS does not have the unique parent element for condition with id " << rCondition.Id() << ".\n";

        const auto& r_parent_element = rCondition.GetValue(NEIGHBOUR_ELEMENTS)[0];
        const auto& r_properties = r_parent_element.GetProperties();

        KRATOS_ERROR_IF_NOT(r_properties.Has(DENSITY)) << "DENSITY is not present in the properties of the parent element with id " << r_parent_element.Id() << " for the condition with id " << rCondition.Id() << ".\n";
        KRATOS_ERROR_IF_NOT(r_properties[DENSITY] > 0.0) << "DENSITY is not greater than zero in the properties of the parent element with id " << r_parent_element.Id() << " for the condition with id " << rCondition.Id() << ".\n";

        KRATOS_ERROR_IF_NOT(r_properties.Has(DYNAMIC_VISCOSITY)) << "DYNAMIC_VISCOSITY is not present in the properties of the parent element with id " << r_parent_element.Id() << " for the condition with id " << rCondition.Id() << ".\n";
        KRATOS_ERROR_IF_NOT(r_properties[DYNAMIC_VISCOSITY] > 0.0) << "DYNAMIC_VISCOSITY is not greater than zero in the properties of the parent element with id " << r_parent_element.Id() << " for the condition with id " << rCondition.Id() << ".\n";
    });

    return 0;

    KRATOS_CATCH("");
}

std::string ComputeYPlusProcess::Info() const
{
    return std::string("ComputeYPlusProcess");
}

void ComputeYPlusProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void ComputeYPlusProcess::PrintData(std::ostream& rOStream) const
{
}

void ComputeYPlusProcess::CalculateConditionReactions()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    VariableUtils().SetNonHistoricalVariableToZero(REACTION, r_model_part.Conditions());

    block_for_each(r_model_part.Conditions(), [](ConditionType& rCondition) {
        auto& r_reaction = rCondition.GetValue(REACTION);

        const auto& r_geometry = rCondition.GetGeometry();
        const double area_fraction = r_geometry.Area() / r_geometry.PointsNumber();

        for (const auto& r_node : r_geometry) {
            const double nodal_area = norm_2(r_node.FastGetSolutionStepValue(NORMAL));
            noalias(r_reaction) += area_fraction * r_node.FastGetSolutionStepValue(REACTION) / nodal_area;
        }
    });

    KRATOS_CATCH("");
}

} // namespace Kratos.
