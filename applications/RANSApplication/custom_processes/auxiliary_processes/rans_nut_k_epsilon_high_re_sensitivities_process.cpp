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

// Application includes
#include "custom_utilities/rans_check_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_nut_k_epsilon_high_re_sensitivities_process.h"

namespace Kratos
{
RansNutKEpsilonHighReSensitivitiesProcess::RansNutKEpsilonHighReSensitivitiesProcess(
    Model& rModel, Parameters& rParameters)
    : mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"      : 0,
            "c_mu"            : 0.09
        })");

    mrParameters.ValidateAndAssignDefaults(default_parameters);

    mEchoLevel = mrParameters["echo_level"].GetInt();
    mModelPartName = mrParameters["model_part_name"].GetString();
    mCmu = mrParameters["c_mu"].GetDouble();

    KRATOS_CATCH("");
}

int RansNutKEpsilonHighReSensitivitiesProcess::Check()
{
    KRATOS_TRY

    RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);

    const ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, TURBULENT_KINETIC_ENERGY);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(
        r_model_part, TURBULENT_ENERGY_DISSIPATION_RATE);

    return 0;

    KRATOS_CATCH("");
}

void RansNutKEpsilonHighReSensitivitiesProcess::ExecuteInitializeSolutionStep()
{
    Execute();
}

void RansNutKEpsilonHighReSensitivitiesProcess::Execute()
{
    KRATOS_TRY

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    ModelPart::NodesContainerType& r_nodes = r_model_part.Nodes();
    int number_of_nodes = r_nodes.size();

#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(r_nodes.begin() + i_node);
        Vector nut_partial_derivatives(2);

        const double& tke = r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
        const double& epsilon =
            r_node.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);

        nut_partial_derivatives[0] = 2.0 * mCmu * tke / epsilon;
        nut_partial_derivatives[1] = -1.0 * mCmu * std::pow(tke / epsilon, 2);

        r_node.SetValue(RANS_NUT_SCALAR_PARTIAL_DERIVATIVES, nut_partial_derivatives);
    }

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 1)
        << "Calculated k-epsilon high Re nu_t sensitivities for nodes in"
        << mModelPartName << "\n";

    KRATOS_CATCH("");
}

std::string RansNutKEpsilonHighReSensitivitiesProcess::Info() const
{
    return std::string("RansNutKEpsilonHighReSensitivitiesProcess");
}

void RansNutKEpsilonHighReSensitivitiesProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansNutKEpsilonHighReSensitivitiesProcess::PrintData(std::ostream& rOStream) const
{
}

} // namespace Kratos.
