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

// Application includes
#include "includes/checks.h"
#include "includes/define.h"

// Include base h
#include "rans_check_utilities.h"

namespace Kratos
{
namespace RansCheckUtilities
{
bool CheckIfModelPartExists(const Model& rModel, const std::string& rModelPartName)
{
    KRATOS_TRY

    if (!rModel.HasModelPart(rModelPartName))
    {
        const std::vector<std::string>& r_model_part_names = rModel.GetModelPartNames();

        std::string msg;
        msg = rModel.Info() + " doesn't have " + rModelPartName +
              ". Available model parts are: \n";
        for (std::string model_part_name : r_model_part_names)
            msg += "     " + model_part_name + "\n";

        KRATOS_ERROR << msg;
    }

    return true;

    KRATOS_CATCH("");
}

template <typename TVariableType>
bool CheckIfVariableExistsInModelPart(const ModelPart& rModelPart, const TVariableType& rVariable)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(!rModelPart.HasNodalSolutionStepVariable(rVariable))
        << rModelPart.Name() << " doesn't have "
        << rVariable.Name() << " in NodalSolutionStepDataContainer. Please add it as a SolutionStepVariable.";

    return true;

    KRATOS_CATCH("");
}

// template instantiations
template bool CheckIfVariableExistsInModelPart(const ModelPart&, const Variable<int>&);
template bool CheckIfVariableExistsInModelPart(const ModelPart&, const Variable<double>&);
template bool CheckIfVariableExistsInModelPart(const ModelPart&,
                                               const Variable<array_1d<double, 3>>&);
} // namespace RansCheckUtilities

} // namespace Kratos
