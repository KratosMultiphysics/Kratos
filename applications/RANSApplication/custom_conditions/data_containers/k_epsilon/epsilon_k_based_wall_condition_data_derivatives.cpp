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
#include "custom_conditions/data_containers/k_epsilon/condition_data_utilities.h"
#include "custom_utilities/fluid_calculation_utilities.h"
#include "custom_utilities/rans_adjoint_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "custom_conditions/data_containers/k_epsilon/epsilon_k_based_wall_condition_data_derivatives.h"

namespace Kratos
{
namespace KEpsilonWallConditionData
{
/***************************************************************/
/********************** Static Operations **********************/
/***************************************************************/

template <unsigned int TDim>
void EpsilonKBasedWallConditionDataDerivatives<TDim>::Data::Check(
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
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA))
        << "TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA is not found in process info.\n";
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
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);

        KRATOS_CHECK_DOF_IN_NODE(RANS_SCALAR_2_ADJOINT_1, r_node);
    }

    KRATOS_CATCH("");
}

/***************************************************************/
/************************ Condition Data ***********************/
/***************************************************************/

template <unsigned int TDim>
EpsilonKBasedWallConditionDataDerivatives<TDim>::Data::Data(
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
    mBeta = r_properties[WALL_SMOOTHNESS_BETA];

    noalias(mUnitNormal) = rGeometry.GetValue(NORMAL);
    mUnitNormal /= norm_2(mUnitNormal);
    mNormalMagnitude = rGeometry.Area();

    // get wall height from the condition data container
    mWallHeight = rGeometry.GetValue(DISTANCE);

    KRATOS_CATCH("");
}

template <unsigned int TDim>
const Variable<double>& EpsilonKBasedWallConditionDataDerivatives<TDim>::Data::GetAdjointScalarVariable()
{
    return RANS_SCALAR_2_ADJOINT_1;
}

template <unsigned int TDim>
void EpsilonKBasedWallConditionDataDerivatives<TDim>::Data::CalculateGaussPointData(
    const Vector& rN,
    const int Step)
{
    FluidCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), rN,
        std::tie(mTurbulentKineticEnergy, TURBULENT_KINETIC_ENERGY),
        std::tie(mWallVelocity, VELOCITY));

    mUTau = mCmu25 * std::sqrt(std::max(mTurbulentKineticEnergy, 0.0));

    mWallVelocityMagnitude = norm_2(mWallVelocity);

    mWallFlux = KEpsilonConditionDataUtilities::CalculateWallFlux(
        mrParameters.mKinematicViscosity, mEpsilonSigma, mUTau,
        mrParameters.mKappa, mrParameters.mYPlus);
}

/***************************************************************/
/********************* Velocity Derivative *********************/
/***************************************************************/

template <unsigned int TDim>
double EpsilonKBasedWallConditionDataDerivatives<TDim>::UDerivative::CalculateWallFluxDerivative(
    const IndexType ConditionNodeIndex,
    const IndexType DirectionIndex,
    const double W,
    const Vector& rN,
    const double WDerivative,
    const double DetJDerivative,
    const IndexType ParentElementNodeIndex) const
{
    double y_plus_derivative = 0.0;
    if (mrData.mWallVelocityMagnitude > 1e-16) {

        const double wall_velocity_magnitude_derivative = mrData.mWallVelocity[DirectionIndex] * rN[ConditionNodeIndex] / mrData.mWallVelocityMagnitude;

        if (mrData.mrParameters.mYPlus > mrData.mYPlusLimit) {
            double temp;
            RansAdjointUtilities::CalculateYPlusAndUtauDerivative(
                y_plus_derivative, temp, mrData.mrParameters.mYPlus,
                mrData.mWallVelocityMagnitude, wall_velocity_magnitude_derivative, mrData.mWallHeight, 0.0, mrData.mrParameters.mKinematicViscosity, mrData.mrParameters.mKappa, mrData.mBeta, mrData.mYPlusLimit);
        }
    }

    return KEpsilonConditionDataUtilities::CalculateWallFluxDerivative(
        mrData.mrParameters.mKinematicViscosity, mrData.mEpsilonSigma, mrData.mUTau, 0.0,
        mrData.mrParameters.mKappa, mrData.mrParameters.mYPlus, y_plus_derivative);
}

template <unsigned int TDim>
double EpsilonKBasedWallConditionDataDerivatives<TDim>::UDerivative::CalculateWallFluxDerivative(
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
double EpsilonKBasedWallConditionDataDerivatives<TDim>::KDerivative::CalculateWallFluxDerivative(
    const IndexType ConditionNodeIndex,
    const IndexType DirectionIndex,
    const double W,
    const Vector& rN,
    const double WDerivative,
    const double DetJDerivative,
    const IndexType ParentElementNodeIndex) const
{
    double u_tau_derivative = 0.0;

    if (mrData.mTurbulentKineticEnergy > 0.0) {
        u_tau_derivative = mrData.mCmu25 * 0.5 * rN[ConditionNodeIndex] /
                           std::sqrt(mrData.mTurbulentKineticEnergy);
    }

    return KEpsilonConditionDataUtilities::CalculateWallFluxDerivative(
        mrData.mrParameters.mKinematicViscosity, mrData.mEpsilonSigma, mrData.mUTau, u_tau_derivative,
        mrData.mrParameters.mKappa, mrData.mrParameters.mYPlus, 0.0);
}

template <unsigned int TDim>
double EpsilonKBasedWallConditionDataDerivatives<TDim>::KDerivative::CalculateWallFluxDerivative(
    const IndexType DirectionIndex,
    const double W,
    const Vector& rN,
    const IndexType ParentElementNodeIndex) const
{
    // zero because wall velocity is computed only with the condition geometry nodes
    return 0.0;
}

/***************************************************************/
/********************** Epsilon Derivative *********************/
/***************************************************************/

template <unsigned int TDim>
double EpsilonKBasedWallConditionDataDerivatives<TDim>::EpsilonDerivative::CalculateWallFluxDerivative(
    const IndexType ConditionNodeIndex,
    const IndexType DirectionIndex,
    const double W,
    const Vector& rN,
    const double WDerivative,
    const double DetJDerivative,
    const IndexType ParentElementNodeIndex) const
{
    return 0.0;
}

template <unsigned int TDim>
double EpsilonKBasedWallConditionDataDerivatives<TDim>::EpsilonDerivative::CalculateWallFluxDerivative(
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
double EpsilonKBasedWallConditionDataDerivatives<TDim>::ShapeDerivative::CalculateWallFluxDerivative(
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
        }
    }

    return KEpsilonConditionDataUtilities::CalculateWallFluxDerivative(
        mrData.mrParameters.mKinematicViscosity, mrData.mEpsilonSigma, mrData.mUTau, 0.0,
        mrData.mrParameters.mKappa, mrData.mrParameters.mYPlus, y_plus_derivative);
}

template <unsigned int TDim>
double EpsilonKBasedWallConditionDataDerivatives<TDim>::ShapeDerivative::CalculateWallFluxDerivative(
    const IndexType DirectionIndex,
    const double W,
    const Vector& rN,
    const IndexType ParentElementNodeIndex) const
{
    // here we compute derivatives w.r.t. parent element nodes
    double y_plus_derivative = 0.0;
    if (mrData.mWallVelocityMagnitude > 1e-16) {
        if (mrData.mrParameters.mYPlus > mrData.mYPlusLimit) {
            const double wall_height_derivative =
                RansAdjointUtilities::CalculateWallHeightParentElementDerivative(
                    mrData.GetGeometry(), mrData.mrParentElementGeometry,
                    DirectionIndex, mrData.mUnitNormal, ZeroVector(3));

            double u_tau_derivative;
            RansAdjointUtilities::CalculateYPlusAndUtauDerivative(
                y_plus_derivative, u_tau_derivative, mrData.mrParameters.mYPlus,
                mrData.mWallVelocityMagnitude, 0.0, mrData.mWallHeight, wall_height_derivative,
                mrData.mrParameters.mKinematicViscosity, mrData.mrParameters.mKappa, mrData.mBeta, mrData.mYPlusLimit);

            // this does not have any effect since mNumberOfGaussPoints is 1.
            // otherwise this may result in undesired behaviour
        }
    }

    return KEpsilonConditionDataUtilities::CalculateWallFluxDerivative(
        mrData.mrParameters.mKinematicViscosity, mrData.mEpsilonSigma, mrData.mUTau, 0.0,
        mrData.mrParameters.mKappa, mrData.mrParameters.mYPlus, y_plus_derivative);
}

// template instantiations

template class EpsilonKBasedWallConditionDataDerivatives<2>;
template class EpsilonKBasedWallConditionDataDerivatives<3>;

} // namespace KEpsilonWallConditionData

} // namespace Kratos