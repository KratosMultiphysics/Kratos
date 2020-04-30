//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Dharmin Shah (https://github.com/sdharmin)
//                   Bence Rochlitz (https://github.com/bencerochlitz)
//
//  Supervised by:   Jordi Cotela (https://github.com/jcotela)
//                   Suneth Warnakulasuriya (https://github.com/sunethwarna)

// System includes
#include <cmath>

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_utilities/rans_check_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_k_omega_wall_function_process.h"
namespace Kratos
{
/// Constructor
RansKOmegaWallFunctionProcess::RansKOmegaWallFunctionProcess(Model& rModel, Parameters& rParameters)
    : mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
    {
        "model_part_name"         : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "turbulent_mixing_length" : 0.005,
        "c_mu"                    : 0.09,
        "echo_level"              : 0,
        "is_fixed"                : true
    })");

    mrParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    mTurbulentMixingLength = mrParameters["turbulent_mixing_length"].GetDouble();
    mIsConstrained = mrParameters["is_fixed"].GetBool();
    mEchoLevel = mrParameters["echo_level"].GetInt();
    mModelPartName = mrParameters["model_part_name"].GetString();
    mCmu_75 = std::pow(mrParameters["c_mu"].GetDouble(), 0.75);

    KRATOS_ERROR_IF(mTurbulentMixingLength < std::numeric_limits<double>::epsilon())
        << "turbulent_mixing_length should be greater than zero.\n";

    KRATOS_CATCH("");
}
RansKOmegaWallFunctionProcess::RansKOmegaWallFunctionProcess(Model& rModel,
                                                             const std::string& rModelPartName)
    : mrModel(rModel), mModelPartName(rModelPartName)
{
}

/// Destructor.
RansKOmegaWallFunctionProcess::~RansKOmegaWallFunctionProcess()
{
}

void RansKOmegaWallFunctionProcess::ExecuteInitialize()
{
    if (mIsConstrained)
    {
        ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

        const int number_of_nodes = r_model_part.NumberOfNodes();
#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_model_part.NodesBegin() + i_node);
            r_node.Fix(TURBULENT_KINETIC_ENERGY);
            r_node.Fix(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
        }

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << "Fixed k_omega dofs in " << mModelPartName << ".\n";
    }

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);
    VariableUtils().SetNonHistoricalVariableToZero(
        NUMBER_OF_NEIGHBOUR_CONDITIONS, r_model_part.Nodes());

    const int number_of_conditions = r_model_part.NumberOfConditions();
#pragma omp parallel for
    for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond)
    {
        ModelPart::ConditionType& r_cond = *(r_model_part.ConditionsBegin() + i_cond);
        ModelPart::ConditionType::GeometryType& r_geometry = r_cond.GetGeometry();
        for (IndexType i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node)
        {
            NodeType& r_node = r_geometry[i_node];
            r_node.SetLock();
            r_node.GetValue(NUMBER_OF_NEIGHBOUR_CONDITIONS) += 1;
            r_node.UnSetLock();
        }
    }

    r_model_part.GetCommunicator().AssembleNonHistoricalData(NUMBER_OF_NEIGHBOUR_CONDITIONS);

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Calculated number of neighbour conditions in " << mModelPartName << ".\n";
}

void RansKOmegaWallFunctionProcess::ExecuteInitializeSolutionStep()
{
    Execute();
}

void RansKOmegaWallFunctionProcess::Execute()
{
    KRATOS_TRY

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    ModelPart::NodesContainerType& r_nodes = r_model_part.Nodes();
    VariableUtils().SetHistoricalVariableToZero(TURBULENT_KINETIC_ENERGY, r_nodes);
    VariableUtils().SetHistoricalVariableToZero(
        TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_nodes);

    ModelPart::ConditionsContainerType& r_conditions = r_model_part.Conditions();
    const int number_of_conditions = r_conditions.size();
    const double mYPlusLimit = 11.06;
    const double mYPlusLowerLimit = 2.0;
    const double mBetaZero = 0.0808;
    const double mCmu25 = std::pow(0.09, 0.25);
    const double mC1 = 0.1;
    const double mKappa = 0.41;

#pragma omp parallel for
    for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond)
    {
        ModelPart::ConditionType& r_condition = *(r_conditions.begin() + i_cond);
        const double y_plus = std::max(r_condition.GetValue(RANS_Y_PLUS), mYPlusLowerLimit);
        const double u_tau = norm_2(r_condition.GetValue(FRICTION_VELOCITY));

        const double omega_vis = 6.0 * std::pow(u_tau / y_plus, 2) / (mBetaZero);
        const double omega_log = std::pow(u_tau / mCmu25, 2) / (mKappa * y_plus);
        const double omega = std::sqrt(std::pow(omega_vis, 2) + std::pow(omega_log, 2));

        const double tke_low = mC1 * std::pow(y_plus * u_tau, 2);
        const double tke_high = std::pow(u_tau / mCmu25, 2);

        ModelPart::ConditionType::GeometryType& r_geometry = r_condition.GetGeometry();
        for (IndexType i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node)
        {
            ModelPart::NodeType& r_node = r_geometry[i_node];
            const double nu = r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);

            r_node.SetLock();
            r_node.FastGetSolutionStepValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE) +=
                omega / nu;
            if (y_plus > mYPlusLimit)
            {
                r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY) += tke_high;
            }
            else
            {
                r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY) += tke_low;
            }
            r_node.UnSetLock();
        }
    }

    r_model_part.GetCommunicator().AssembleCurrentData(TURBULENT_KINETIC_ENERGY);
    r_model_part.GetCommunicator().AssembleCurrentData(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);

    const int number_of_nodes = r_nodes.size();
#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        ModelPart::NodeType& r_node = *(r_nodes.begin() + i_node);
        const double count = r_node.GetValue(NUMBER_OF_NEIGHBOUR_CONDITIONS);

        double& omega = r_node.FastGetSolutionStepValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
        omega = std::max(omega / count, 1e-12);
        double& tke = r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
        tke = std::max(tke / count, 1e-12);
    }

    KRATOS_CATCH("");
}

int RansKOmegaWallFunctionProcess::Check()
{
    RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, TURBULENT_KINETIC_ENERGY);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(
        r_model_part, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);

    return 0;
}

/// Turn back information as a string.
std::string RansKOmegaWallFunctionProcess::Info() const
{
    return std::string("RansKOmegaWallFunctionProcess");
}

/// Print information about this object.
void RansKOmegaWallFunctionProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

/// Print object's data.
void RansKOmegaWallFunctionProcess::PrintData(std::ostream& rOStream) const
{
}

void RansKOmegaWallFunctionProcess::CalculateTurbulentValues(NodeType& rNode)
{
    const double tke = rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
    rNode.FastGetSolutionStepValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE) =
        std::sqrt(std::max(tke, 0.0)) / (mCmu_75 * mTurbulentMixingLength);
}

} // namespace Kratos.