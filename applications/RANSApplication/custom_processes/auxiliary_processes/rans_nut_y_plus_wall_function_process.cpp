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
#include "utilities/variable_utils.h"

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

    return 0;

    KRATOS_CATCH("");
}

void RansNutYPlusWallFunctionProcess::ExecuteInitialize()
{
    CalculateConditionNeighbourCount();
}

void RansNutYPlusWallFunctionProcess::CalculateConditionNeighbourCount()
{
    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);
    VariableUtils().SetNonHistoricalVariableToZero(
        NUMBER_OF_NEIGHBOUR_CONDITIONS, r_model_part.Nodes());

    const int number_of_conditions = r_model_part.NumberOfConditions();
#pragma omp parallel for
    for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond)
    {
        ConditionType& r_cond = *(r_model_part.ConditionsBegin() + i_cond);
        if (RansCalculationUtilities::IsWall(r_cond))
        {
            ConditionGeometryType& r_geometry = r_cond.GetGeometry();
            for (IndexType i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node)
            {
                NodeType& r_node = r_geometry[i_node];
                r_node.SetLock();
                r_node.GetValue(NUMBER_OF_NEIGHBOUR_CONDITIONS) += 1;
                r_node.UnSetLock();
            }
        }
    }

    r_model_part.GetCommunicator().AssembleNonHistoricalData(NUMBER_OF_NEIGHBOUR_CONDITIONS);
}

void RansNutYPlusWallFunctionProcess::Execute()
{
    KRATOS_TRY

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);
    const int number_of_nodes = r_model_part.NumberOfNodes();
#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(r_model_part.NodesBegin() + i_node);
        if (RansCalculationUtilities::IsWall(r_node))
        {
            r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY) = 0.0;
        }
    }

    const int number_of_conditions = r_model_part.NumberOfConditions();
#pragma omp parallel for
    for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond)
    {
        ConditionType& r_cond = *(r_model_part.ConditionsBegin() + i_cond);
        if (RansCalculationUtilities::IsWall(r_cond))
        {
            ConditionType::GeometryType& r_geometry = r_cond.GetGeometry();
            const double y_plus = r_cond.GetValue(RANS_Y_PLUS);

            for (IndexType i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node)
            {
                NodeType& r_node = r_geometry[i_node];
                const double nu = r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);
                r_node.SetLock();
                r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY) +=
                    mVonKarman * y_plus * nu;
                r_node.UnSetLock();
            }
        }
    }
    r_model_part.GetCommunicator().AssembleCurrentData(TURBULENT_VISCOSITY);

#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(r_model_part.NodesBegin() + i_node);
        if (RansCalculationUtilities::IsWall(r_node))
        {
            double& r_nut = r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY);
            const double number_of_neighbour_conditions =
                static_cast<double>(r_node.GetValue(NUMBER_OF_NEIGHBOUR_CONDITIONS));
            r_nut = RansCalculationUtilities::SoftMax(
                r_nut / number_of_neighbour_conditions, mMinValue);
            r_node.FastGetSolutionStepValue(VISCOSITY) =
                r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY) + r_nut;
        }
    }

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
