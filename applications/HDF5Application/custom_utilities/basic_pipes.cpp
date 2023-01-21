//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Máté Kelemen
//

// Internal includes
#include "basic_pipes.h"


namespace Kratos::Pipes {


ModelPartFromModel::ModelPartFromModel(std::string&& rModelPartName) noexcept
    : mModelPartName(std::move(rModelPartName))
{
}


ModelPartFromModel::ModelPartFromModel(const std::string& rModelPartName)
    : mModelPartName(rModelPartName)
{
}


ModelPartFromModel::ModelPartFromModel(const Parameters& rParameters)
    : mModelPartName()
{
    KRATOS_TRY
    KRATOS_ERROR_IF_NOT(rParameters.Has("model_part_name"))
        << "ModelPartFromModel requires a \"model_part_name\" entry in the input parameters but found none in:\n"
        << rParameters;

    const auto model_part_name = rParameters["model_part_name"];
    KRATOS_ERROR_IF_NOT(model_part_name.IsString()) << "ModelPartFromModel expects \"model_part_name\" as a string, but got\n"
        << model_part_name;

    mModelPartName = model_part_name.GetString();
    KRATOS_CATCH("");
}


TimeFromProcessInfo::TimeFromProcessInfo()
    : VariableFromProcessInfo<decltype(TIME)>(TIME)
{
}


TimeFromProcessInfo::TimeFromProcessInfo(const Parameters& rParameters)
    : TimeFromProcessInfo()
{
}


StepFromProcessInfo::StepFromProcessInfo()
    : VariableFromProcessInfo<decltype(STEP)>(STEP)
{
}


StepFromProcessInfo::StepFromProcessInfo(const Parameters& rParameters)
    : StepFromProcessInfo()
{
}


} // namespace Kratos::Pipes
