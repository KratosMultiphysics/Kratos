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

// External includes

// Project includes
#include "includes/define.h"

// Application includes

// Include base h
#include "rans_vtk_output_process.h"
namespace Kratos
{
/// Constructor
RansVTKOutputProcess::RansVTKOutputProcess(Model& rModel, Parameters& rParameters)
    : mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mEchoLevel = rParameters["echo_level"].GetInt();
    mModelPartName = rParameters["model_part_name"].GetString();

    if (rParameters["output_settings"].Has("model_part_name")) {
        KRATOS_WARNING_IF(this->Info(),
                          rParameters["output_settings"]["model_part_name"].GetString() != mModelPartName)
            << "VTKOutput is specified with "
            << rParameters["output_settings"]["model_part_name"].GetString()
            << ", and " << mModelPartName << " is specified as main model part. Using "
            << mModelPartName << " for VTKOutput.\n";
    } else {
        rParameters["output_settings"].AddEmptyValue("model_part_name");
    }
    rParameters["output_settings"]["model_part_name"].SetString(mModelPartName);
    mIndex = 0;

    mpVTKOutput = Kratos::make_shared<VtkOutput>(
        rModel.GetModelPart(mModelPartName), rParameters["output_settings"]);
    mOutputPath = rParameters["output_settings"]["output_path"].GetString();

    this->UpdateExecutionPointsList(rParameters["execution_points"].GetStringArray());

    KRATOS_CATCH("");
}

void RansVTKOutputProcess::ExecuteOperation()
{
    KRATOS_TRY

    std::stringstream file_name;

    file_name << mModelPartName << "_output_itr_" << mIndex++;

    mpVTKOutput->PrintOutput(file_name.str());

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Written VTK output to " << mOutputPath << "/" << file_name.str() << ".\n";

    KRATOS_CATCH("");
}

/// Turn back information as a string.
std::string RansVTKOutputProcess::Info() const
{
    return std::string("RansVTKOutputProcess");
}

/// Print information about this object.
void RansVTKOutputProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

/// Print object's data.
void RansVTKOutputProcess::PrintData(std::ostream& rOStream) const
{
}

const Parameters RansVTKOutputProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
    {
        "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "echo_level"      : 0,
        "execution_points": [],
        "output_settings" : {}
    })");

    return default_parameters;
}

} // namespace Kratos.