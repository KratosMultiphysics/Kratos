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
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_check_utilities.h"
#include "rans_modelling_application_variables.h"

// Include base h
#include "rans_nut_k_epsilon_high_re_calculation_process.h"

namespace Kratos
{
RansNutKEpsilonHighReCalculationProcess::RansNutKEpsilonHighReCalculationProcess(Model& rModel, Parameters rParameters)
    : mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"      : 0,
            "c_mu"            : 0.09,
            "min_nut_factor"  : 1e-4
        })");

    mrParameters.ValidateAndAssignDefaults(default_parameters);

    mEchoLevel = mrParameters["echo_level"].GetInt();
    mModelPartName = mrParameters["model_part_name"].GetString();
    mCmu = mrParameters["c_mu"].GetDouble();
    mMinNutFactor = mrParameters["min_nut_factor"].GetDouble();

    KRATOS_ERROR_IF(mMinNutFactor < 0.0)
        << "\"min_nut_factor\"=" << mMinNutFactor
        << " is negative. Please provide a positive value.\n";

    KRATOS_CATCH("");
}

RansNutKEpsilonHighReCalculationProcess::~RansNutKEpsilonHighReCalculationProcess()
{
}

int RansNutKEpsilonHighReCalculationProcess::Check()
{
    KRATOS_TRY

    RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);

    const ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, TURBULENT_KINETIC_ENERGY);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(
        r_model_part, TURBULENT_ENERGY_DISSIPATION_RATE);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, TURBULENT_VISCOSITY);

    return 0;

    KRATOS_CATCH("");
}

void RansNutKEpsilonHighReCalculationProcess::Execute()
{
    KRATOS_TRY

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    NodesContainerType& r_nodes = r_model_part.Nodes();
    int number_of_nodes = r_nodes.size();

    // #pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(r_nodes.begin() + i_node);
        const double tke = r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
        const double epsilon =
            r_node.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
        const double nu_t = EvmKepsilonModelUtilities::CalculateTurbulentViscosity(
            mCmu, tke, epsilon, 1.0);
        const double nu = r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);
        const double min_nut_value = nu * mMinNutFactor;
        const double soft_nu_t = RansCalculationUtilities::SoftMax(nu_t, min_nut_value);

        // KRATOS_ERROR_IF(soft_nu_t == 0.0)
        //     << "TURBULENT_VISCOSITY=" << nu_t << "  at node with id=" << r_node.Id()
        //     << ", located at=" << r_node.Coordinates()
        //     << " is less than the machine precision [ TURBULENT_VISCOSITY < "
        //     << std::numeric_limits<double>::epsilon() << " ] with soft minimum bound of min_nu_t="
        //     << min_nut_value << ". This leads to numerical instabilities in the proceeding calculations. "
        //     << "Please try increasing \"min_nut_factor\" to increase min_nu_t. "
        //        "[ min_nu_t = min_nut_factor * KINEMATIC_VISCOSITY, "
        //        "KINEMATIC_VISCOSITY="
        //     << nu << ", min_nu_t_factor=" << mMinNutFactor << " ].\n";

        r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY) = soft_nu_t;
    }

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 1)
        << "Calculated nu_t for nodes in" << mModelPartName << "\n";

    KRATOS_CATCH("");
}

std::string RansNutKEpsilonHighReCalculationProcess::Info() const
{
    return std::string("RansNutKEpsilonHighReCalculationProcess");
}

void RansNutKEpsilonHighReCalculationProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansNutKEpsilonHighReCalculationProcess::PrintData(std::ostream& rOStream) const
{
}

} // namespace Kratos.
