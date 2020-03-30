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

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/define.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_check_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_nut_y_plus_wall_function_process.h"

namespace Kratos
{
RansNutYPlusWallFunctionProcess::RansNutYPlusWallFunctionProcess(Model& rModel, Parameters rParameters)
    : mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"      : 0,
            "c_mu"            : 0.09,
            "von_karman"      : 0.41,
            "beta"            : 5.2,
            "min_value"       : 1e-18
        })");

    mrParameters.ValidateAndAssignDefaults(default_parameters);

    mEchoLevel = mrParameters["echo_level"].GetInt();
    mModelPartName = mrParameters["model_part_name"].GetString();
    mCmu = mrParameters["c_mu"].GetDouble();
    mVonKarman = mrParameters["von_karman"].GetDouble();
    mBeta = mrParameters["beta"].GetDouble();
    mMinValue = mrParameters["min_value"].GetDouble();
    mLimitYPlus =
        RansCalculationUtilities::CalculateLogarithmicYPlusLimit(mVonKarman, mBeta);

    KRATOS_CATCH("");
}

int RansNutYPlusWallFunctionProcess::Check()
{
    KRATOS_TRY

    RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);

    const ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, KINEMATIC_VISCOSITY);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, TURBULENT_VISCOSITY);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, RANS_Y_PLUS);

    return 0;

    KRATOS_CATCH("");
}

void RansNutYPlusWallFunctionProcess::Execute()
{
    KRATOS_TRY

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    const int number_of_nodes = r_model_part.NumberOfNodes();

    unsigned int number_of_modified_nu_t_nodes = 0;

#pragma omp parallel for reduction(+ : number_of_modified_nu_t_nodes)
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(r_model_part.NodesBegin() + i_node);
        const double y_plus = r_node.FastGetSolutionStepValue(RANS_Y_PLUS);
        const double nu = r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);

        if (y_plus > mLimitYPlus)
        {
            r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY) = mVonKarman * y_plus * nu;
            ++number_of_modified_nu_t_nodes;
        }
        else
        {
            r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY) = mMinValue;
        }
    }

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Applied nu_t y_plus wall function to "
        << number_of_modified_nu_t_nodes << " of total "
        << r_model_part.NumberOfNodes() << " nodes in " << mModelPartName << "\n";

    KRATOS_CATCH("");
}

std::string RansNutYPlusWallFunctionProcess::Info() const
{
    return std::string("RansNutYPlusWallFunctionProcess");
}

void RansNutYPlusWallFunctionProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansNutYPlusWallFunctionProcess::PrintData(std::ostream& rOStream) const
{
}

} // namespace Kratos.
