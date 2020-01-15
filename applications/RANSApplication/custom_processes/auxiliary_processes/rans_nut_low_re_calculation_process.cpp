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
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_utilities.h"
#include "custom_utilities/rans_check_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_nut_low_re_calculation_process.h"

namespace Kratos
{
RansNutLowReCalculationProcess::RansNutLowReCalculationProcess(Model& rModel, Parameters rParameters)
    : mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"      : 0,
            "c_mu"            : 0.09,
            "min_value"       : 1e-15
        })");

    mrParameters.ValidateAndAssignDefaults(default_parameters);

    mEchoLevel = mrParameters["echo_level"].GetInt();
    mModelPartName = mrParameters["model_part_name"].GetString();
    mCmu = mrParameters["c_mu"].GetDouble();
    mMinValue = mrParameters["min_value"].GetDouble();

    KRATOS_CATCH("");
}

int RansNutLowReCalculationProcess::Check()
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

void RansNutLowReCalculationProcess::Execute()
{
    KRATOS_TRY

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    const ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
    const double c_mu = r_process_info[TURBULENCE_RANS_C_MU];

    NodesContainerType& r_nodes = r_model_part.Nodes();
    int number_of_nodes = r_nodes.size();

#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(r_nodes.begin() + i_node);
        const double epsilon =
            r_node.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
        if (epsilon > 0.0)
        {
            const double tke = r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
            const double y_plus = r_node.FastGetSolutionStepValue(RANS_Y_PLUS);
            const double f_mu = EvmKepsilonModelUtilities::CalculateFmu(y_plus);
            const double nu_t = EvmKepsilonModelUtilities::CalculateTurbulentViscosity(
                c_mu, tke, epsilon, f_mu);

            r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY) = nu_t;
        }
        else
        {
            r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY) = mMinValue;
        }
    }

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 1)
        << "Calculated nu_t for nodes in" << mModelPartName << "\n";

    KRATOS_CATCH("");
}

std::string RansNutLowReCalculationProcess::Info() const
{
    return std::string("RansNutLowReCalculationProcess");
}

void RansNutLowReCalculationProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansNutLowReCalculationProcess::PrintData(std::ostream& rOStream) const
{
}

} // namespace Kratos.
