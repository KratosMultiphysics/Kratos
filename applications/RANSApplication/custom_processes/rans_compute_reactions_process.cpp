//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/constitutive_law.h"
#include "includes/define.h"
#include "includes/node.h"
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_utilities/fluid_calculation_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"

// Include base h
#include "rans_compute_reactions_process.h"

namespace Kratos
{

RansComputeReactionsProcess::RansComputeReactionsProcess(
    Model& rModel,
    Parameters rParameters)
: mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mModelPartName = rParameters["model_part_name"].GetString();
    mEchoLevel = rParameters["echo_level"].GetInt();
    this->UpdateExecutionPointsList(rParameters["execution_points"].GetStringArray());

    KRATOS_CATCH("");
}

RansComputeReactionsProcess::RansComputeReactionsProcess(
    Model& rModel,
    const std::string& rModelPartName,
    const std::vector<std::string>& rExecutionPoints,
    const int EchoLevel)
    : mrModel(rModel),
      mModelPartName(rModelPartName),
      mEchoLevel(EchoLevel)
{
    this->UpdateExecutionPointsList(rExecutionPoints);
}

int RansComputeReactionsProcess::Check()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    block_for_each(r_model_part.Conditions(), [](ConditionType& rCondition){
        KRATOS_ERROR_IF_NOT(rCondition.Has(NEIGHBOUR_ELEMENTS))
            << "NEIGHBOUR_ELEMENTS is not found in the condition data container. [ Condition.Id() = "
            << rCondition.Id() << " ].\n";

        KRATOS_ERROR_IF(rCondition.GetValue(NEIGHBOUR_ELEMENTS).size() != 1)
            << "NEIGHBOUR_ELEMENTS in condition contains its parent, which needs to be only one. [ NEIGHBOUR_ELEMENTS.size() = "
            << rCondition.GetValue(NEIGHBOUR_ELEMENTS).size() << ", Condition.Id() = " << rCondition.Id() << " ].\n";

        KRATOS_ERROR_IF_NOT(rCondition.Has(NORMAL))
            << "NORMAL is not found in the condition data container. [ Condition.Id() = "
            << rCondition.Id() << " ].\n";
    });

    return 0;

    KRATOS_CATCH("");
}

void RansComputeReactionsProcess::Initialize()
{
    KRATOS_TRY

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0) << "Initialized compute reactions process on " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

std::string RansComputeReactionsProcess::Info() const
{
    return std::string("RansComputeReactionsProcess");
}

void RansComputeReactionsProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansComputeReactionsProcess::PrintData(std::ostream& rOStream) const
{
}

const Parameters RansComputeReactionsProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "execution_points": ["after_coupling_solve_step"],
            "echo_level"      : 0
        })");

    return default_parameters;
}

void RansComputeReactionsProcess::ExecuteOperation()
{
    KRATOS_TRY

    if (!mIsInitialized) {
        Initialize();
    }

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);
    const auto& r_process_info = r_model_part.GetProcessInfo();

    if (r_process_info[DOMAIN_SIZE] == 2) {
        CalculateReactions<2>(r_model_part);
    } else if (r_process_info[DOMAIN_SIZE] == 3) {
        CalculateReactions<3>(r_model_part);
    } else {
        KRATOS_ERROR << "Unsupported DOMAIN_SIZE [ DOMAIN_SIZE " << r_process_info[DOMAIN_SIZE] << " ].\n";
    }

    VariableUtils().SetHistoricalVariableToZero(REACTION, r_model_part.Nodes());

    block_for_each(r_model_part.Conditions(), [](ConditionType& rCondition){
        const int number_of_nodes = rCondition.GetGeometry().PointsNumber();
        const auto& r_reaction = rCondition.GetValue(REACTION);
        for (auto& r_node : rCondition.GetGeometry()) {
            r_node.SetLock();
            r_node.FastGetSolutionStepValue(REACTION) += r_reaction / number_of_nodes;
            r_node.UnSetLock();
        }
    });

    r_model_part.GetCommunicator().AssembleCurrentData(REACTION);

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0) << "Calculated reactions on " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

template<>
void RansComputeReactionsProcess::CalculateStrainRate<2>(
    Vector& rStrainRate,
    const GeometryType& rElementGeometry,
    const Matrix& rdNdX) const
{
    noalias(rStrainRate) = ZeroVector(3);
    for (IndexType i = 0; i < rElementGeometry.PointsNumber(); ++i) {
        const auto& r_velocity = rElementGeometry[i].FastGetSolutionStepValue(VELOCITY);
        rStrainRate[0] += rdNdX(i, 0) * r_velocity[0];
        rStrainRate[1] += rdNdX(i, 1) * r_velocity[1];
        rStrainRate[2] += rdNdX(i, 0) * r_velocity[1] + rdNdX(i, 1) * r_velocity[0];
    }
}

template<>
void RansComputeReactionsProcess::CalculateStrainRate<3>(
    Vector& rStrainRate,
    const GeometryType& rElementGeometry,
    const Matrix& rdNdX) const
{
    noalias(rStrainRate) = ZeroVector(6);
    for (IndexType i = 0; i < rElementGeometry.PointsNumber(); ++i) {
        const auto& r_velocity = rElementGeometry[i].FastGetSolutionStepValue(VELOCITY);
        rStrainRate[0] += rdNdX(i, 0) * r_velocity[0];
        rStrainRate[1] += rdNdX(i, 1) * r_velocity[1];
        rStrainRate[2] += rdNdX(i, 2) * r_velocity[2];
        rStrainRate[3] += rdNdX(i, 0) * r_velocity[1] + rdNdX(i, 1) * r_velocity[0];
        rStrainRate[4] += rdNdX(i, 1) * r_velocity[2] + rdNdX(i, 2) * r_velocity[1];
        rStrainRate[5] += rdNdX(i, 0) * r_velocity[2] + rdNdX(i, 2) * r_velocity[0];
    }
}

template<>
void RansComputeReactionsProcess::CalculateViscousStressTensorReactionContribution<2>(
    array_1d<double, 3>& rReaction,
    const Vector& rViscousStress,
    const array_1d<double, 3>& rNormal) const
{
    // viscous stress is defined in the following order
    // sigma_xx, sigma_yy, sigma_xy

    rReaction[0] = rViscousStress[0] * rNormal[0] + rViscousStress[2] * rNormal[1];
    rReaction[1] = rViscousStress[2] * rNormal[0] + rViscousStress[1] * rNormal[1];
    rReaction[2] = 0.0;
}

template<>
void RansComputeReactionsProcess::CalculateViscousStressTensorReactionContribution<3>(
    array_1d<double, 3>& rReaction,
    const Vector& rViscousStress,
    const array_1d<double, 3>& rNormal) const
{
    // viscous stress is defined in the following order
    // sigma_xx(0), sigma_yy(1), sigma_zz(2), sigma_xy(3), sigma_yz(4), sigma_xz(5)

    rReaction[0] = rViscousStress[0] * rNormal[0] + rViscousStress[3] * rNormal[1] + rViscousStress[5] * rNormal[2];
    rReaction[1] = rViscousStress[3] * rNormal[0] + rViscousStress[1] * rNormal[1] + rViscousStress[4] * rNormal[2];
    rReaction[2] = rViscousStress[5] * rNormal[0] + rViscousStress[4] * rNormal[1] + rViscousStress[2] * rNormal[2];
}

template<unsigned int TDim>
void RansComputeReactionsProcess::CalculateReactions(ModelPart& rModelPart) const
{
    KRATOS_TRY

    constexpr IndexType strain_size = (TDim - 1) * 3;

    struct TLSType
    {
        Vector mParentWs;
        Matrix mParentNs;
        ShapeFunctionDerivativesArrayType mParentdNdXs;

        Vector mConditionWs;
        Matrix mConditionNs;

        Vector mStrainRate = ZeroVector(strain_size);
        Vector mViscousStress = ZeroVector(strain_size);
    };

    const auto& r_process_info = rModelPart.GetProcessInfo();

    block_for_each(rModelPart.Conditions(), TLSType(), [&](ConditionType& rCondition, TLSType& rTLS) {
        // first calculate elemental data to get the stress tensor
        const auto& r_parent_element = rCondition.GetValue(NEIGHBOUR_ELEMENTS)[0];
        const auto& r_parent_properties = r_parent_element.GetProperties();
        const auto& r_parent_geometry = r_parent_element.GetGeometry();

        RansCalculationUtilities::CalculateGeometryData(
            r_parent_geometry,
            GeometryData::IntegrationMethod::GI_GAUSS_1, rTLS.mParentWs, rTLS.mParentNs, rTLS.mParentdNdXs);

        const Vector& parent_N = row(rTLS.mParentNs, 0);
        const Matrix& parent_dNdX = rTLS.mParentdNdXs[0];

        // calculate stress tensor
        auto constitutive_law = r_parent_element.GetProperties().GetValue(CONSTITUTIVE_LAW);

        auto cl_values = ConstitutiveLaw::Parameters(
            r_parent_geometry, r_parent_properties, r_process_info);

        Flags& cl_options = cl_values.GetOptions();
        cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS);

        cl_values.SetShapeFunctionsValues(parent_N);
        cl_values.SetShapeFunctionsDerivatives(parent_dNdX);

        CalculateStrainRate<TDim>(rTLS.mStrainRate, r_parent_geometry, parent_dNdX);

        cl_values.SetStrainVector(rTLS.mStrainRate);   //this is the input parameter
        cl_values.SetStressVector(rTLS.mViscousStress);  //this is an ouput parameter

        constitutive_law->CalculateMaterialResponseCauchy(cl_values);

        // now calculate on the condition
        const auto& r_condition_geometry = rCondition.GetGeometry();

        RansCalculationUtilities::CalculateConditionGeometryData(
            r_condition_geometry, GeometryData::IntegrationMethod::GI_GAUSS_1, rTLS.mConditionWs, rTLS.mConditionNs);

        // get the normal (this includes the area of the condition as well)
        auto r_normal = r_condition_geometry.GetValue(NORMAL);
        r_normal *= r_condition_geometry.Area() / norm_2(r_normal);

        array_1d<double, 3> reaction;

        // calculate viscous stress reaction contribution
        CalculateViscousStressTensorReactionContribution<TDim>(
            reaction, rTLS.mViscousStress, r_normal);

        // calculate pressure contribution
        double pressure;
        FluidCalculationUtilities::EvaluateInPoint(
            r_condition_geometry, row(rTLS.mConditionNs, 0), std::tie(pressure, PRESSURE));
        noalias(reaction) -= r_normal * pressure;

        // now set the elemental reaction
        rCondition.GetValue(REACTION) = reaction;
    });

    KRATOS_CATCH("");
}

template void RansComputeReactionsProcess::CalculateReactions<2>(ModelPart&) const;
template void RansComputeReactionsProcess::CalculateReactions<3>(ModelPart&) const;

} // namespace Kratos.
