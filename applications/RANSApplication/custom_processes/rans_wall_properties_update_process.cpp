//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Dharmin Shah
//                   Bence Rochlitz
//
//  Supervised by:   Jordi Cotela
//                   Suneth Warnakulasuriya
//

// System includes
#include <cmath>

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"
#include "utilities/normal_calculation_utils.h"

// Application includes
#include "custom_utilities/fluid_calculation_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_wall_properties_update_process.h"
namespace Kratos
{
/// Constructor
RansWallPropertiesUpdateProcess::RansWallPropertiesUpdateProcess(
    Model& rModel,
    Parameters& rParameters)
    : mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mEchoLevel = rParameters["echo_level"].GetInt();
    mModelPartName = rParameters["model_part_name"].GetString();
    mUpdateWallNormals = rParameters["update_wall_normals"].GetBool();
    mUpdateConditionWallHeights = rParameters["update_condition_wall_heights"].GetBool();
    this->UpdateExecutionPointsList(rParameters["update_execution_points"].GetStringArray());

    KRATOS_CATCH("");
}

RansWallPropertiesUpdateProcess::RansWallPropertiesUpdateProcess(
    Model& rModel,
    const std::string& rModelPartName,
    const bool UpdateWallNormals,
    const bool UpdateConditionWallHeights,
    const std::vector<std::string>& rUpdateExecutionPoints,
    const int EchoLevel)
    : mrModel(rModel),
      mModelPartName(rModelPartName),
      mEchoLevel(EchoLevel),
      mUpdateWallNormals(UpdateWallNormals),
      mUpdateConditionWallHeights(UpdateConditionWallHeights)
{
    KRATOS_TRY

    this->UpdateExecutionPointsList(rUpdateExecutionPoints);

    KRATOS_CATCH("");
}

int RansWallPropertiesUpdateProcess::Check()
{
    KRATOS_TRY

    return 0;

    KRATOS_CATCH("");
}

/// Turn back information as a string.
std::string RansWallPropertiesUpdateProcess::Info() const
{
    return std::string("RansWallPropertiesUpdateProcess");
}

/// Print information about this object.
void RansWallPropertiesUpdateProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

/// Print object's data.
void RansWallPropertiesUpdateProcess::PrintData(std::ostream& rOStream) const
{
}

const Parameters RansWallPropertiesUpdateProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
    {
        "model_part_name"              : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "echo_level"                   : 0,
        "update_execution_points"      : ["initialize"],
        "update_condition_wall_heights": true,
        "update_wall_normals"          : true
    })");

    return default_parameters;
}

void RansWallPropertiesUpdateProcess::ExecuteOperation()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    if (mUpdateWallNormals) {
        NormalCalculationUtils().CalculateOnSimplex(r_model_part);
        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0) << "Updated normals in " << mModelPartName << ".\n";
    }

    if (mUpdateConditionWallHeights) {
        block_for_each(r_model_part.Conditions(), [&](ModelPart::ConditionType& rCondition) {
            const array_1d<double, 3>& r_normal = rCondition.GetValue(NORMAL);

            const double wall_height = RansCalculationUtilities::CalculateWallHeight(rCondition, r_normal);

            KRATOS_ERROR_IF(wall_height == 0.0) << this->Info() << " has zero wall height.\n";

            rCondition.SetValue(DISTANCE, wall_height);
        });

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0) << "Updated condition wall heights in " << mModelPartName << ".\n";
    }

    KRATOS_CATCH("");
}

} // namespace Kratos.