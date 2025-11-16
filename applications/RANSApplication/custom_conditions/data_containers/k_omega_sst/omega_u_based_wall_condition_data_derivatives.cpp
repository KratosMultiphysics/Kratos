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
#include "includes/checks.h"
#include "containers/variable.h"
#include "geometries/geometry.h"
#include "includes/cfd_variables.h"
#include "includes/node.h"
#include "includes/properties.h"
#include "includes/ublas_interface.h"

// Application includes
#include "custom_conditions/data_containers/k_omega/condition_data_utilities.h"
#include "custom_elements/data_containers/k_omega_sst/element_data_derivative_utilities.h"
#include "custom_elements/data_containers/k_omega_sst/element_data_utilities.h"
#include "custom_utilities/fluid_calculation_utilities.h"
#include "custom_utilities/rans_adjoint_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "omega_u_based_wall_condition_data_derivatives.h"

namespace Kratos
{
namespace KOmegaSSTWallConditionData
{
/***************************************************************/
/********************** Static Operations **********************/
/***************************************************************/

template <unsigned int TDim>
void OmegaUBasedWallConditionDataDerivatives<TDim>::Data::Check(
    const Condition& rCondition,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = rCondition.GetGeometry();
    const auto& r_condition_properties = rCondition.GetProperties();
    const int number_of_nodes = r_geometry.PointsNumber();

    KRATOS_ERROR_IF_NOT(r_geometry.GetValue(RANS_IS_WALL_FUNCTION_ACTIVE) == 1)
        << "Wall function check is called for a condition with non-activated wall function. [ Condition.Id() = "
        << rCondition.Id() << " ].\n";

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENCE_RANS_C_MU))
        << "TURBULENCE_RANS_C_MU is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_1))
        << "TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_1 is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_2))
        << "TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_2 is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(WALL_CORRECTION_FACTOR))
        << "WALL_CORRECTION_FACTOR is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(VON_KARMAN))
        << "VON_KARMAN is not found in process info.\n";

    KRATOS_ERROR_IF_NOT(r_condition_properties.Has(RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT))
        << "RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT is not found in condition properties [ Condition.Id() = "
        << rCondition.Id() << ", Properties.Id() = " << r_condition_properties.Id() << " ].\n";
    KRATOS_ERROR_IF_NOT(r_condition_properties.Has(WALL_SMOOTHNESS_BETA))
        << "WALL_SMOOTHNESS_BETA is not found in condition properties [ Condition.Id() = "
        << rCondition.Id() << ", Properties.Id() = " << r_condition_properties.Id() << " ].\n";

    KRATOS_ERROR_IF_NOT(norm_2(r_geometry.GetValue(NORMAL)) > 0.0)
        << "NORMAL is not found or not properly initialized in condition [ Condition.Id() = " << rCondition.Id()
        << " ].\n";
    KRATOS_ERROR_IF_NOT(r_geometry.Has(NORMAL_SHAPE_DERIVATIVE))
        << "NORMAL_SHAPE_DERIVATIVE is not found in condition [ Condition.Id() = " << rCondition.Id()
        << " ].\n";

    KRATOS_ERROR_IF_NOT(r_geometry.GetValue(DISTANCE) > 0.0)
    << "DISTANCE is not found or not properly initialized in condition [ Condition.Id() = " << rCondition.Id()
    << " ].\n";

    for (int i_node = 0; i_node < number_of_nodes; ++i_node) {
        const auto& r_node = r_geometry[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);

        KRATOS_CHECK_DOF_IN_NODE(RANS_SCALAR_2_ADJOINT_1, r_node);
    }

    KRATOS_CATCH("");
}

/***************************************************************/
/************************ Condition Data ***********************/
/***************************************************************/

template <unsigned int TDim>
OmegaUBasedWallConditionDataDerivatives<TDim>::Data::Data(
    const GeometryType& rGeometry,
    const Properties& rProperties,
    const ProcessInfo& rProcessInfo,
    const ScalarWallFluxConditionData::Parameters& rParameters)
    : BaseType(rGeometry, rProperties, rProcessInfo),
      mrParentElementGeometry(rGeometry.GetValue(NEIGHBOUR_ELEMENTS)[0].GetGeometry()),
      mrParameters(rParameters)
{
    KRATOS_TRY

    CalculateConstants(rProcessInfo);

    const auto& r_properties = this->GetConditionProperties();
    mYPlusLimit = r_properties[RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT];

    noalias(mUnitNormal) = rGeometry.GetValue(NORMAL);
    mUnitNormal /= norm_2(mUnitNormal);
    mNormalMagnitude = rGeometry.Area();

    // get wall height from the condition data container
    mWallHeight = rGeometry.GetValue(DISTANCE);

    KRATOS_CATCH("");
}

template <unsigned int TDim>
const Variable<double>& OmegaUBasedWallConditionDataDerivatives<TDim>::Data::GetAdjointScalarVariable()
{
    return RANS_SCALAR_2_ADJOINT_1;
}

template <unsigned int TDim>
void OmegaUBasedWallConditionDataDerivatives<TDim>::Data::CalculateGaussPointData(
    const Vector& rN,
    const int Step)
{
    FluidCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), rN,
        std::tie(mTurbulentKineticEnergy, TURBULENT_KINETIC_ENERGY),
        std::tie(mTurbulentSpecificEnergyDissipationRate, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE),
        std::tie(mWallVelocity, VELOCITY));

    mWallVelocityMagnitude = norm_2(mWallVelocity);

    mUTau = mWallVelocityMagnitude / ((1.0 / mrParameters.mKappa) * std::log(mrParameters.mYPlus) + mBeta);

    mF1 = KOmegaSSTElementData::CalculateF1(
        mTurbulentKineticEnergy, mTurbulentSpecificEnergyDissipationRate, mrParameters.mKinematicViscosity, mWallHeight, mBetaStar, 1e-10, mSigmaOmega2);

    mBlendedSigma = KOmegaSSTElementData::CalculateBlendedPhi(mSigmaOmega1, mSigmaOmega2, mF1) * mWallCorrectionFactor;

    mWallFlux = KOmegaConditionDataUtilities::CalculateWallFlux(
        mrParameters.mKinematicViscosity, mBlendedSigma, mUTau, mCmu25,
        mrParameters.mKappa, mrParameters.mYPlus);
}

/***************************************************************/
/********************* Velocity Derivative *********************/
/***************************************************************/

template <unsigned int TDim>
double OmegaUBasedWallConditionDataDerivatives<TDim>::UDerivative::CalculateWallFluxDerivative(
    const IndexType ConditionNodeIndex,
    const IndexType DirectionIndex,
    const double W,
    const Vector& rN,
    const double WDerivative,
    const double DetJDerivative,
    const IndexType ParentElementNodeIndex) const
{
    double y_plus_derivative = 0.0;
    double u_tau_derivative = 0.0;
    if (mrData.mWallVelocityMagnitude > 1e-16) {
        const double wall_velocity_magnitude_derivative = mrData.mWallVelocity[DirectionIndex] * rN[ConditionNodeIndex] / mrData.mWallVelocityMagnitude;

        const double denominator = (1.0 / mrData.mrParameters.mKappa) * std::log(mrData.mrParameters.mYPlus) + mrData.mBeta;
        u_tau_derivative += wall_velocity_magnitude_derivative / denominator;

        if (mrData.mrParameters.mYPlus > mrData.mYPlusLimit) {
            double temp;
            RansAdjointUtilities::CalculateYPlusAndUtauDerivative(
                y_plus_derivative, temp, mrData.mrParameters.mYPlus,
                mrData.mWallVelocityMagnitude, wall_velocity_magnitude_derivative, mrData.mWallHeight, 0.0, mrData.mrParameters.mKinematicViscosity, mrData.mrParameters.mKappa, mrData.mBeta, mrData.mYPlusLimit);

            const double denominator_derivative = (1.0 / mrData.mrParameters.mKappa) * y_plus_derivative / mrData.mrParameters.mYPlus;
            u_tau_derivative -= mrData.mWallVelocityMagnitude * denominator_derivative / std::pow(denominator, 2);
        }
    }

    const double f_1_derivative = KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateF1Derivative(
        mrData.mTurbulentKineticEnergy, 0.0, mrData.mTurbulentSpecificEnergyDissipationRate,
        0.0, mrData.mrParameters.mKinematicViscosity, mrData.mWallHeight, 0.0, mrData.mBetaStar,
        1e-10, 0.0, mrData.mSigmaOmega2);

    const double blended_sigma_derivative =
        KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateBlendedPhiDerivative(
            mrData.mSigmaOmega1, 0.0, mrData.mSigmaOmega2, 0.0, mrData.mF1, f_1_derivative) *
        mrData.mWallCorrectionFactor;

    return KOmegaConditionDataUtilities::CalculateWallFluxDerivative(
        mrData.mrParameters.mKinematicViscosity, mrData.mBlendedSigma, blended_sigma_derivative, mrData.mUTau, u_tau_derivative,
        mrData.mCmu25, mrData.mrParameters.mKappa, mrData.mrParameters.mYPlus, y_plus_derivative);
}

template <unsigned int TDim>
double OmegaUBasedWallConditionDataDerivatives<TDim>::UDerivative::CalculateWallFluxDerivative(
    const IndexType DirectionIndex,
    const double W,
    const Vector& rN,
    const IndexType ParentElementNodeIndex) const
{
    return 0.0;
}

/***************************************************************/
/************************* K Derivative ************************/
/***************************************************************/

template <unsigned int TDim>
double OmegaUBasedWallConditionDataDerivatives<TDim>::KDerivative::CalculateWallFluxDerivative(
    const IndexType ConditionNodeIndex,
    const IndexType DirectionIndex,
    const double W,
    const Vector& rN,
    const double WDerivative,
    const double DetJDerivative,
    const IndexType ParentElementNodeIndex) const
{
    const double f_1_derivative = KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateF1Derivative(
        mrData.mTurbulentKineticEnergy, rN[ConditionNodeIndex], mrData.mTurbulentSpecificEnergyDissipationRate,
        0.0, mrData.mrParameters.mKinematicViscosity, mrData.mWallHeight, 0.0, mrData.mBetaStar,
        1e-10, 0.0, mrData.mSigmaOmega2);

    const double blended_sigma_derivative =
        KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateBlendedPhiDerivative(
            mrData.mSigmaOmega1, 0.0, mrData.mSigmaOmega2, 0.0, mrData.mF1, f_1_derivative) *
        mrData.mWallCorrectionFactor;

    return KOmegaConditionDataUtilities::CalculateWallFluxDerivative(
        mrData.mrParameters.mKinematicViscosity, mrData.mBlendedSigma, blended_sigma_derivative, mrData.mUTau, 0.0,
        mrData.mCmu25, mrData.mrParameters.mKappa, mrData.mrParameters.mYPlus, 0.0);
}

template <unsigned int TDim>
double OmegaUBasedWallConditionDataDerivatives<TDim>::KDerivative::CalculateWallFluxDerivative(
    const IndexType DirectionIndex,
    const double W,
    const Vector& rN,
    const IndexType ParentElementNodeIndex) const
{
    // zero because wall velocity is computed only with the condition geometry nodes
    return 0.0;
}

/***************************************************************/
/********************** Omega Derivative *********************/
/***************************************************************/

template <unsigned int TDim>
double OmegaUBasedWallConditionDataDerivatives<TDim>::OmegaDerivative::CalculateWallFluxDerivative(
    const IndexType ConditionNodeIndex,
    const IndexType DirectionIndex,
    const double W,
    const Vector& rN,
    const double WDerivative,
    const double DetJDerivative,
    const IndexType ParentElementNodeIndex) const
{
    const double f_1_derivative = KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateF1Derivative(
        mrData.mTurbulentKineticEnergy, 0.0, mrData.mTurbulentSpecificEnergyDissipationRate,
        rN[ConditionNodeIndex], mrData.mrParameters.mKinematicViscosity, mrData.mWallHeight, 0.0, mrData.mBetaStar,
        1e-10, 0.0, mrData.mSigmaOmega2);

    const double blended_sigma_derivative =
        KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateBlendedPhiDerivative(
            mrData.mSigmaOmega1, 0.0, mrData.mSigmaOmega2, 0.0, mrData.mF1, f_1_derivative) *
        mrData.mWallCorrectionFactor;

    return KOmegaConditionDataUtilities::CalculateWallFluxDerivative(
        mrData.mrParameters.mKinematicViscosity, mrData.mBlendedSigma, blended_sigma_derivative, mrData.mUTau, 0.0,
        mrData.mCmu25, mrData.mrParameters.mKappa, mrData.mrParameters.mYPlus, 0.0);
}

template <unsigned int TDim>
double OmegaUBasedWallConditionDataDerivatives<TDim>::OmegaDerivative::CalculateWallFluxDerivative(
    const IndexType DirectionIndex,
    const double W,
    const Vector& rN,
    const IndexType ParentElementNodeIndex) const
{
    // zero because wall velocity is computed only with the condition geometry nodes
    return 0.0;
}

/***************************************************************/
/*********************** Shape Derivative **********************/
/***************************************************************/

template <unsigned int TDim>
double OmegaUBasedWallConditionDataDerivatives<TDim>::ShapeDerivative::CalculateWallFluxDerivative(
    const IndexType ConditionNodeIndex,
    const IndexType DirectionIndex,
    const double W,
    const Vector& rN,
    const double WDerivative,
    const double DetJDerivative,
    const IndexType ParentElementNodeIndex) const
{
    // here we compute derivatives w.r.t. condition nodes
    double y_plus_derivative = 0.0;
    double u_tau_derivative = 0.0;
    double blended_sigma_derivative = 0.0;
    if (mrData.mWallVelocityMagnitude > 1e-16) {
        if (mrData.mrParameters.mYPlus > mrData.mYPlusLimit) {
            const Vector& temp = row(mrData.GetGeometry().GetValue(NORMAL_SHAPE_DERIVATIVE), ConditionNodeIndex * TDim + DirectionIndex);
            array_1d<double, 3> normal_derivative = ZeroVector(3);
            for (IndexType i = 0; i < TDim; ++i) {
                normal_derivative[i] = temp[i];
            }

            const auto& unit_normal_derivative = RansAdjointUtilities::CalculateUnitVectorDerivative(
                mrData.mNormalMagnitude, mrData.mUnitNormal, normal_derivative);

            const double wall_height_derivative =
                RansAdjointUtilities::CalculateWallHeightConditionDerivative(
                    mrData.GetGeometry(), mrData.mrParentElementGeometry,
                    DirectionIndex, mrData.mUnitNormal, unit_normal_derivative);

            double value;
            RansAdjointUtilities::CalculateYPlusAndUtauDerivative(
                y_plus_derivative, value, mrData.mrParameters.mYPlus,
                mrData.mWallVelocityMagnitude, 0.0, mrData.mWallHeight, wall_height_derivative,
                mrData.mrParameters.mKinematicViscosity, mrData.mrParameters.mKappa, mrData.mBeta, mrData.mYPlusLimit);

            const double f_1_derivative = KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateF1Derivative(
                mrData.mTurbulentKineticEnergy, 0.0, mrData.mTurbulentSpecificEnergyDissipationRate,
                0.0, mrData.mrParameters.mKinematicViscosity, mrData.mWallHeight, wall_height_derivative, mrData.mBetaStar,
                1e-10, 0.0, mrData.mSigmaOmega2);

            blended_sigma_derivative =
                    KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateBlendedPhiDerivative(
                        mrData.mSigmaOmega1, 0.0, mrData.mSigmaOmega2, 0.0, mrData.mF1, f_1_derivative) *
                    mrData.mWallCorrectionFactor;

            const double denominator = (1.0 / mrData.mrParameters.mKappa) * std::log(mrData.mrParameters.mYPlus) + mrData.mBeta;
            const double denominator_derivative = (1.0 / mrData.mrParameters.mKappa) * y_plus_derivative / mrData.mrParameters.mYPlus;
            u_tau_derivative -= mrData.mWallVelocityMagnitude * denominator_derivative / std::pow(denominator, 2);
        }
    }

    return KOmegaConditionDataUtilities::CalculateWallFluxDerivative(
        mrData.mrParameters.mKinematicViscosity, mrData.mBlendedSigma, blended_sigma_derivative, mrData.mUTau, u_tau_derivative,
        mrData.mCmu25, mrData.mrParameters.mKappa, mrData.mrParameters.mYPlus, y_plus_derivative);
}

template <unsigned int TDim>
double OmegaUBasedWallConditionDataDerivatives<TDim>::ShapeDerivative::CalculateWallFluxDerivative(
    const IndexType DirectionIndex,
    const double W,
    const Vector& rN,
    const IndexType ParentElementNodeIndex) const
{
    // here we compute derivatives w.r.t. parent element nodes
    double y_plus_derivative = 0.0;
    double u_tau_derivative = 0.0;
    double blended_sigma_derivative = 0.0;
    if (mrData.mWallVelocityMagnitude > 1e-16) {
        if (mrData.mrParameters.mYPlus > mrData.mYPlusLimit) {
            const double wall_height_derivative =
                RansAdjointUtilities::CalculateWallHeightParentElementDerivative(
                    mrData.GetGeometry(), mrData.mrParentElementGeometry,
                    DirectionIndex, mrData.mUnitNormal, ZeroVector(3));

            double value;
            RansAdjointUtilities::CalculateYPlusAndUtauDerivative(
                y_plus_derivative, value, mrData.mrParameters.mYPlus,
                mrData.mWallVelocityMagnitude, 0.0, mrData.mWallHeight, wall_height_derivative,
                mrData.mrParameters.mKinematicViscosity, mrData.mrParameters.mKappa, mrData.mBeta, mrData.mYPlusLimit);

            const double f_1_derivative = KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateF1Derivative(
                mrData.mTurbulentKineticEnergy, 0.0, mrData.mTurbulentSpecificEnergyDissipationRate,
                0.0, mrData.mrParameters.mKinematicViscosity, mrData.mWallHeight, wall_height_derivative, mrData.mBetaStar,
                1e-10, 0.0, mrData.mSigmaOmega2);

            blended_sigma_derivative =
                    KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateBlendedPhiDerivative(
                        mrData.mSigmaOmega1, 0.0, mrData.mSigmaOmega2, 0.0, mrData.mF1, f_1_derivative) *
                    mrData.mWallCorrectionFactor;

            const double denominator = (1.0 / mrData.mrParameters.mKappa) * std::log(mrData.mrParameters.mYPlus) + mrData.mBeta;
            const double denominator_derivative = (1.0 / mrData.mrParameters.mKappa) * y_plus_derivative / mrData.mrParameters.mYPlus;
            u_tau_derivative -= mrData.mWallVelocityMagnitude * denominator_derivative / std::pow(denominator, 2);
        }
    }

    return KOmegaConditionDataUtilities::CalculateWallFluxDerivative(
        mrData.mrParameters.mKinematicViscosity, mrData.mBlendedSigma, blended_sigma_derivative, mrData.mUTau, u_tau_derivative,
        mrData.mCmu25, mrData.mrParameters.mKappa, mrData.mrParameters.mYPlus, y_plus_derivative);
}

// template instantiations

template class OmegaUBasedWallConditionDataDerivatives<2>;
template class OmegaUBasedWallConditionDataDerivatives<3>;

} // namespace KOmegaSSTWallConditionData

} // namespace Kratos