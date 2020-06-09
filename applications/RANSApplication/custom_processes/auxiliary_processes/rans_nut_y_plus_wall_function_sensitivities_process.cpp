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

// System includes
#include <string>

// External includes

// Project includes
#include "custom_utilities/rans_check_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_nut_y_plus_wall_function_sensitivities_process.h"

namespace Kratos
{
RansNutYPlusWallFunctionSensitivitiesProcess::RansNutYPlusWallFunctionSensitivitiesProcess(
    Model& rModel, Parameters& rParameters)
    : mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"      : 0
        })");

    mrParameters.ValidateAndAssignDefaults(default_parameters);

    mEchoLevel = mrParameters["echo_level"].GetInt();
    mModelPartName = mrParameters["model_part_name"].GetString();

    KRATOS_CATCH("");
}

int RansNutYPlusWallFunctionSensitivitiesProcess::Check()
{
    KRATOS_TRY

    RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);

    return 0;

    KRATOS_CATCH("");
}

void RansNutYPlusWallFunctionSensitivitiesProcess::ExecuteInitializeSolutionStep()
{
    Execute();
}

void RansNutYPlusWallFunctionSensitivitiesProcess::Execute()
{
    KRATOS_TRY

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    const int number_of_nodes = r_model_part.NumberOfNodes();

    unsigned int number_of_modified_nu_t_nodes = 0;

#pragma omp parallel for reduction(+ : number_of_modified_nu_t_nodes)
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(r_model_part.NodesBegin() + i_node);
        Vector nut_partial_derivatives(2);
        nut_partial_derivatives.clear();
        r_node.SetValue(RANS_NUT_SCALAR_PARTIAL_DERIVATIVES, nut_partial_derivatives);
        number_of_modified_nu_t_nodes++;
    }

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Applied nu_t y_plus wall function sensitivities to "
        << number_of_modified_nu_t_nodes << " of total "
        << r_model_part.NumberOfNodes() << " nodes in " << mModelPartName << "\n";

    KRATOS_CATCH("");
}

std::string RansNutYPlusWallFunctionSensitivitiesProcess::Info() const
{
    return std::string("RansNutYPlusWallFunctionSensitivitiesProcess");
}

void RansNutYPlusWallFunctionSensitivitiesProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansNutYPlusWallFunctionSensitivitiesProcess::PrintData(std::ostream& rOStream) const
{
}

} // namespace Kratos.
