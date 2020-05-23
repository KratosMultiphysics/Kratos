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
#include "includes/cfd_variables.h"

// Application includes
#include "custom_utilities/rans_check_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_logarithmic_y_plus_velocity_sensitivities_process.h"

namespace Kratos
{
RansLogarithmicYPlusVelocitySensitivitiesProcess::RansLogarithmicYPlusVelocitySensitivitiesProcess(
    Model& rModel, Parameters& rParameters)
    : mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"      : 0,
            "constants": {
                "von_karman"  : 0.41,
                "beta"        : 5.2
            }
        })");

    mrParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    mModelPartName = mrParameters["model_part_name"].GetString();
    mEchoLevel = mrParameters["echo_level"].GetInt();

    mVonKarman = mrParameters["constants"]["von_karman"].GetDouble();
    mBeta = mrParameters["constants"]["beta"].GetDouble();

    KRATOS_CATCH("");
}

int RansLogarithmicYPlusVelocitySensitivitiesProcess::Check()
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

void RansLogarithmicYPlusVelocitySensitivitiesProcess::ExecuteInitializeSolutionStep()
{
    this->Execute();
}

void RansLogarithmicYPlusVelocitySensitivitiesProcess::Execute()
{
    KRATOS_TRY

    ModelPart& model_part = mrModel.GetModelPart(mModelPartName);

    ModelPart::ElementsContainerType& r_elements = model_part.Elements();

    const int number_of_elements = r_elements.size();

    const int domain_size = model_part.GetProcessInfo()[DOMAIN_SIZE];

    const double inv_kappa = 1.0 / mVonKarman;

#pragma omp parallel for
    for (int i_element = 0; i_element < number_of_elements; ++i_element)
    {
        ElementType& r_element = *(r_elements.begin() + i_element);
        GeometryType& r_geometry = r_element.GetGeometry();
        const int number_of_nodes = r_geometry.PointsNumber();

        Matrix r_adjoint_y_plus_matrix(number_of_nodes, domain_size);
        r_adjoint_y_plus_matrix.clear();

        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            const NodeType& r_node = r_geometry[i_node];
            const double y_plus = r_node.FastGetSolutionStepValue(RANS_Y_PLUS);
            const double wall_distance = r_node.FastGetSolutionStepValue(DISTANCE);
            const array_1d<double, 3> velocity = r_node.FastGetSolutionStepValue(VELOCITY);
            const double velocity_magnitude = norm_2(velocity);
            const double nu = r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);

            double value = 0.0;

            if (y_plus > 11.06)
            {
                value = (inv_kappa * std::log(y_plus) + mBeta);
                value = value / (std::pow(value, 2) + velocity_magnitude * wall_distance *
                                                          inv_kappa / (nu * y_plus));
            }
            else
            {
                value = 1.0 / (2.0 * y_plus);
            }
            if (velocity_magnitude > std::numeric_limits<double>::epsilon() &&
                y_plus > std::numeric_limits<double>::epsilon())
            {
                for (int i_dim = 0; i_dim < domain_size; ++i_dim)
                {
                    r_adjoint_y_plus_matrix(i_node, i_dim) =
                        (wall_distance / nu) * value * velocity[i_dim] / velocity_magnitude;
                }
            }
        }
        r_element.SetValue(RANS_Y_PLUS_VELOCITY_DERIVATIVES, r_adjoint_y_plus_matrix);
    }

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "RANS_Y_PLUS_VELOCITY_DERIVATIVES calculated for "
        << r_elements.size() << " elements in " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

std::string RansLogarithmicYPlusVelocitySensitivitiesProcess::Info() const
{
    return std::string("RansLogarithmicYPlusVelocitySensitivitiesProcess");
}

void RansLogarithmicYPlusVelocitySensitivitiesProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansLogarithmicYPlusVelocitySensitivitiesProcess::PrintData(std::ostream& rOStream) const
{
}
} // namespace Kratos
