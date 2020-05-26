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
#include "includes/define.h"

// Application includes
#include "custom_utilities/rans_check_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_epsilon_turbulent_mixing_inlet_process.h"

namespace Kratos
{
RansEpsilonTurbulentMixingLengthInletProcess::RansEpsilonTurbulentMixingLengthInletProcess(
    Model& rModel, Parameters rParameters)
    : mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name"         : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "turbulent_mixing_length" : 0.005,
            "c_mu"                    : 0.09,
            "echo_level"              : 0,
            "is_fixed"                : true,
            "min_value"               : 1e-14
        })");

    mrParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    mTurbulentMixingLength = mrParameters["turbulent_mixing_length"].GetDouble();
    mIsConstrained = mrParameters["is_fixed"].GetBool();
    mEchoLevel = mrParameters["echo_level"].GetInt();
    mModelPartName = mrParameters["model_part_name"].GetString();
    mCmu_75 = std::pow(mrParameters["c_mu"].GetDouble(), 0.75);
    mMinValue = mrParameters["min_value"].GetDouble();

    KRATOS_ERROR_IF(mTurbulentMixingLength < std::numeric_limits<double>::epsilon())
        << "turbulent_mixing_length should be greater than zero.\n";

    KRATOS_ERROR_IF(mMinValue < 0.0) << "Minimum turbulent energy dissipation "
                                        "rate needs to be positive in the "
                                        "modelpart "
                                     << mModelPartName << "\n.";

    KRATOS_CATCH("");
}

void RansEpsilonTurbulentMixingLengthInletProcess::ExecuteInitialize()
{
    if (mIsConstrained)
    {
        ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

        const int number_of_nodes = r_model_part.NumberOfNodes();
#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_model_part.NodesBegin() + i_node);
            r_node.Fix(TURBULENT_ENERGY_DISSIPATION_RATE);
        }

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << "Fixed TURBULENT_ENERGY_DISSIPATION_RATE dofs in "
            << mModelPartName << ".\n";
    }
}

void RansEpsilonTurbulentMixingLengthInletProcess::ExecuteInitializeSolutionStep()
{
    Execute();
}

void RansEpsilonTurbulentMixingLengthInletProcess::Execute()
{
    KRATOS_TRY

    ModelPart::NodesContainerType& r_nodes =
        mrModel.GetModelPart(mModelPartName).Nodes();
    const int number_of_nodes = r_nodes.size();

#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(r_nodes.begin() + i_node);
        CalculateTurbulentValues(r_node);
    }

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Applied epsilon values to " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

int RansEpsilonTurbulentMixingLengthInletProcess::Check()
{
    RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, TURBULENT_KINETIC_ENERGY);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(
        r_model_part, TURBULENT_ENERGY_DISSIPATION_RATE);

    return 0;
}

std::string RansEpsilonTurbulentMixingLengthInletProcess::Info() const
{
    return std::string("RansEpsilonTurbulentMixingLengthInletProcess");
}

void RansEpsilonTurbulentMixingLengthInletProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansEpsilonTurbulentMixingLengthInletProcess::PrintData(std::ostream& rOStream) const
{
}

void RansEpsilonTurbulentMixingLengthInletProcess::CalculateTurbulentValues(NodeType& rNode)
{
    const double tke = rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
    rNode.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE) = std::max(
        mCmu_75 * std::pow(std::max(tke, 0.0), 1.5) / mTurbulentMixingLength, mMinValue);
}

} // namespace Kratos.
