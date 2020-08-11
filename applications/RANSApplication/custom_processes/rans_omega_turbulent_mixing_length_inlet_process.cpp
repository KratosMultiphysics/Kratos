//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Dharmin Shah
//                   Bence Rochlitz
//
//  Supervised by:   Jordi Cotela
//                   Suneth Warnakulasuriya
//

// System includes
#include <cmath>

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "rans_application_variables.h"

// Include base h
#include "rans_omega_turbulent_mixing_length_inlet_process.h"
namespace Kratos
{
/// Constructor
RansOmegaTurbulentMixingLengthInletProcess::RansOmegaTurbulentMixingLengthInletProcess(
    Model& rModel,
    Parameters& rParameters)
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
        "min_value"               : 1e-12
    })");

    mrParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    mTurbulentMixingLength = mrParameters["turbulent_mixing_length"].GetDouble();
    mIsConstrained = mrParameters["is_fixed"].GetBool();
    mEchoLevel = mrParameters["echo_level"].GetInt();
    mModelPartName = mrParameters["model_part_name"].GetString();
    mCmu_25 = std::pow(mrParameters["c_mu"].GetDouble(), 0.25);
    mMinValue = mrParameters["min_value"].GetDouble();

    KRATOS_ERROR_IF(mTurbulentMixingLength < std::numeric_limits<double>::epsilon())
        << "turbulent_mixing_length should be greater than zero.\n";

    KRATOS_ERROR_IF(mMinValue < 0.0) << "Minimum turbulent energy dissipation "
                                        "rate needs to be positive in the "
                                        "modelpart "
                                     << mModelPartName << "\n.";

    KRATOS_CATCH("");
}
/// Destructor.
RansOmegaTurbulentMixingLengthInletProcess::~RansOmegaTurbulentMixingLengthInletProcess()
{
}

void RansOmegaTurbulentMixingLengthInletProcess::ExecuteInitialize()
{
    if (mIsConstrained) {
        auto& r_model_part = mrModel.GetModelPart(mModelPartName);

        const int number_of_nodes = r_model_part.NumberOfNodes();

        BlockPartition<ModelPart::NodesContainerType>(r_model_part.Nodes())
            .for_each([&](ModelPart::NodeType& rNode) {
                rNode.Fix(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
            });

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << "Fixed TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE dofs in "
            << mModelPartName << ".\n";
    }
}

void RansOmegaTurbulentMixingLengthInletProcess::ExecuteInitializeSolutionStep()
{
    Execute();
}

void RansOmegaTurbulentMixingLengthInletProcess::Execute()
{
    KRATOS_TRY

    auto& r_nodes = mrModel.GetModelPart(mModelPartName).Nodes();

    BlockPartition<ModelPart::NodesContainerType>(r_nodes).for_each([&](ModelPart::NodeType& rNode) {
        const double tke = rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
        rNode.FastGetSolutionStepValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE) =
            std::max(std::sqrt(std::max(tke, 0.0)) / (mCmu_25 * mTurbulentMixingLength),
                     mMinValue);
    });

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Applied omega values to " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

int RansOmegaTurbulentMixingLengthInletProcess::Check()
{
    KRATOS_TRY

    const auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY))
        << "TURBULENT_KINETIC_ENERGY is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE))
        << "TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    return 0;

    KRATOS_CATCH("");
}

/// Turn back information as a string.
std::string RansOmegaTurbulentMixingLengthInletProcess::Info() const
{
    return std::string("RansOmegaTurbulentMixingLengthInletProcess");
}

/// Print information about this object.
void RansOmegaTurbulentMixingLengthInletProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

/// Print object's data.
void RansOmegaTurbulentMixingLengthInletProcess::PrintData(std::ostream& rOStream) const
{
}

} // namespace Kratos.