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
#include <cmath>
#include <limits>

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/define.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_check_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_logarithmic_y_plus_calculation_process.h"

namespace Kratos
{
RansLogarithmicYPlusCalculationProcess::RansLogarithmicYPlusCalculationProcess(Model& rModel,
                                                                               Parameters rParameters)
    : mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"      : 0,
            "step"            : 0,
            "max_iterations"  : 100,
            "tolerance"       : 1e-6,
            "constants": {
                "von_karman"  : 0.41,
                "beta"        : 5.2
            }
        })");

    mrParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    mModelPartName = mrParameters["model_part_name"].GetString();

    mEchoLevel = mrParameters["echo_level"].GetInt();
    mStep = mrParameters["step"].GetInt();

    mMaxIterations = mrParameters["max_iterations"].GetInt();
    mTolerance = mrParameters["tolerance"].GetDouble();

    mVonKarman = mrParameters["constants"]["von_karman"].GetDouble();
    mBeta = mrParameters["constants"]["beta"].GetDouble();

    mLimitYPlus =
        RansCalculationUtilities::CalculateLogarithmicYPlusLimit(mVonKarman, mBeta);

    KRATOS_CATCH("");
}

int RansLogarithmicYPlusCalculationProcess::Check()
{
    KRATOS_TRY


    RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);

    const ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, VELOCITY);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, DISTANCE);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, KINEMATIC_VISCOSITY);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, RANS_Y_PLUS);

    return 0;

    KRATOS_CATCH("");
}

void RansLogarithmicYPlusCalculationProcess::Execute()
{
    ModelPart::NodesContainerType& r_nodes =
        mrModel.GetModelPart(mModelPartName).Nodes();

    const int number_of_nodes = r_nodes.size();

#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(r_nodes.begin() + i_node);
        CalculateLogarithmicWallLawYplus(r_node, VELOCITY);
    }

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "RANS_Y_PLUS calculated for nodes in " << mModelPartName << ".\n";
}

std::string RansLogarithmicYPlusCalculationProcess::Info() const
{
    return std::string("RansLogarithmicYPlusCalculationProcess");
}

void RansLogarithmicYPlusCalculationProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansLogarithmicYPlusCalculationProcess::PrintData(std::ostream& rOStream) const
{
}

void RansLogarithmicYPlusCalculationProcess::CalculateLogarithmicWallLawYplus(
    NodeType& rNode, const Variable<array_1d<double, 3>>& rVelocityVariable)
{
    KRATOS_TRY

    const double kinematic_viscosity =
        rNode.FastGetSolutionStepValue(KINEMATIC_VISCOSITY, mStep);
    const double wall_distance = rNode.FastGetSolutionStepValue(DISTANCE, mStep);
    const array_1d<double, 3>& velocity =
        rNode.FastGetSolutionStepValue(rVelocityVariable, mStep);
    const double velocity_magnitude = norm_2(velocity);

    if (velocity_magnitude < std::numeric_limits<double>::epsilon())
    {
        rNode.FastGetSolutionStepValue(RANS_Y_PLUS) = 0.0;
        return;
    }

    KRATOS_ERROR_IF(wall_distance < std::numeric_limits<double>::epsilon())
        << "DISTANCE at node " << rNode.Coordinates()
        << " with id=" << rNode.Id() << " has zero value. Please specify DISTANCE value > 0.0 for all the nodes in "
        << mModelPartName << " to calculate RANS_Y_PLUS.\n";

    // linear region
    double utau = std::sqrt(velocity_magnitude * kinematic_viscosity / wall_distance);
    double yplus = wall_distance * utau / kinematic_viscosity;
    const double inv_von_karman = 1.0 / mVonKarman;

    // log region
    if (yplus > mLimitYPlus)
    {
        unsigned int iter = 0;
        double dx = 1e10;
        const double tol = mTolerance;
        double uplus = inv_von_karman * std::log(yplus) + mBeta;

        while (iter < mMaxIterations && std::fabs(dx) > tol * utau)
        {
            // Newton-Raphson iteration
            double f = utau * uplus - velocity_magnitude;
            double df = uplus + inv_von_karman;
            dx = f / df;

            // Update variables
            utau -= dx;
            yplus = wall_distance * utau / kinematic_viscosity;
            uplus = inv_von_karman * std::log(yplus) + mBeta;
            ++iter;
        }

        KRATOS_WARNING_IF("RansLogarithmicYPlusCalculationProcess", iter == mMaxIterations && mEchoLevel > 0)
            << "Y plus calculation in Wall (logarithmic region) "
               "Newton-Raphson did not converge. "
               "residual > tolerance [ "
            << std::scientific << dx << " > " << std::scientific << tol << " ]\n";
    }
    rNode.FastGetSolutionStepValue(RANS_Y_PLUS) = yplus;

    KRATOS_CATCH("");
}
} // namespace Kratos.