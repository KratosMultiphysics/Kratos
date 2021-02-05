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
    const auto& r_element = rCondition.GetValue(NEIGHBOUR_ELEMENTS)[0];
    const auto& r_element_properties = r_element.GetProperties();
    const auto& r_condition_properties = rCondition.GetProperties();
    const int number_of_nodes = r_geometry.PointsNumber();

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENCE_RANS_C_MU))
        << "TURBULENCE_RANS_C_MU is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA))
        << "TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(r_element_properties.Has(DYNAMIC_VISCOSITY))
        << "DYNAMIC_VISCOSITY is not found in element properties [ Element.Id() = "
        << r_element.Id() << ", Properties.Id() = " << r_element_properties.Id() << " ].\n";
    KRATOS_ERROR_IF_NOT(r_element_properties.Has(DENSITY))
        << "DENSITY is not found in element properties [ Element.Id() = "
        << r_element.Id() << ", Properties.Id() = " << r_element_properties.Id() << " ].\n";
    KRATOS_ERROR_IF_NOT(r_condition_properties.Has(RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT))
        << "RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT is not found in condition properties [ Condition.Id() = "
        << rCondition.Id() << ", Properties.Id() = " << r_condition_properties.Id() << " ].\n";
    KRATOS_ERROR_IF_NOT(r_condition_properties.Has(WALL_SMOOTHNESS_BETA))
        << "WALL_SMOOTHNESS_BETA is not found in condition properties [ Condition.Id() = "
        << rCondition.Id() << ", Properties.Id() = " << r_condition_properties.Id() << " ].\n";
    KRATOS_ERROR_IF_NOT(r_geometry.Has(NORMAL))
        << "NORMAL is not found in condition [ Condition.Id() = " << rCondition.Id()
        << " ].\n";
    KRATOS_ERROR_IF_NOT(r_geometry.Has(NORMAL_SHAPE_DERIVATIVE))
        << "NORMAL_SHAPE_DERIVATIVE is not found in condition [ Condition.Id() = " << rCondition.Id()
        << " ].\n";

    for (int i_node = 0; i_node < number_of_nodes; ++i_node) {
        const auto& r_node = r_geometry[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);;

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
    const ProcessInfo& rProcessInfo)
    : BaseType(rGeometry, rProperties, rProcessInfo),
      mrParentElementGeometry(rGeometry.GetValue(NEIGHBOUR_ELEMENTS)[0].GetGeometry())
{
    KRATOS_TRY

    CalculateConstants(rProcessInfo);

    mCmu = rProcessInfo[TURBULENCE_RANS_C_MU];

    const auto& r_properties = this->GetConditionProperties();
    mYPlusLimit = r_properties[RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT];
    mBeta = r_properties[WALL_SMOOTHNESS_BETA];

    const auto& r_normal = rGeometry.GetValue(NORMAL);
    mNormalMagnitude = norm_2(r_normal);
    mUnitNormal = r_normal / mNormalMagnitude;

    // calculate wall height
    mWallHeight = inner_prod(rGeometry.Center() - mrParentElementGeometry.Center(), mUnitNormal);

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
    auto& cl_parameters = this->GetConstitutiveLawParameters();
    cl_parameters.SetShapeFunctionsValues(rN);

    this->GetConstitutiveLaw().CalculateValue(cl_parameters, EFFECTIVE_VISCOSITY, mKinematicViscosity);
    mKinematicViscosity /= mDensity;

    FluidCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), rN,
        std::tie(mTurbulentKinematicViscosity, TURBULENT_VISCOSITY),
        std::tie(mTurbulentKineticEnergy, TURBULENT_KINETIC_ENERGY),
        std::tie(mWallVelocity, VELOCITY));

    mUTau = mCmu25 * std::sqrt(std::max(mTurbulentKineticEnergy, 0.0));

    mWallVelocityMagnitude = norm_2(mWallVelocity);

    mWallFlux =
        (mKinematicViscosity + mTurbulentKinematicViscosity / mEpsilonSigma) *
        std::pow(mUTau, 5) / (mKappa * std::pow(mYPlus * mKinematicViscosity, 2));
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

        if (mrData.mYPlus > mrData.mYPlusLimit) {
        double u_tau_derivative;
        RansAdjointUtilities::CalculateYPlusAndUtauDerivative(
            y_plus_derivative, u_tau_derivative, mrData.mYPlus,
            mrData.mWallVelocityMagnitude, wall_velocity_magnitude_derivative, mrData.mWallHeight, 0.0, mrData.mKinematicViscosity, mrData.mKappa, mrData.mBeta, mrData.mYPlusLimit);
        }

        // this does not have any effect since mNumberOfGaussPoints is 1.
        // otherwise this may result in undesired behaviour
    }

    return (mrData.mKinematicViscosity +  mrData.mTurbulentKinematicViscosity / mrData.mEpsilonSigma) * std::pow(mrData.mUTau, 5) *
           (-2.0 * y_plus_derivative / std::pow(mrData.mYPlus, 3)) / (mrData.mKappa * std::pow(mrData.mKinematicViscosity, 2));
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
    double nu_t_derivative = 0.0;

    if (mrData.mTurbulentKineticEnergy > 0.0) {
        u_tau_derivative = mrData.mCmu25 * 0.5 * rN[ConditionNodeIndex] /
                           std::sqrt(mrData.mTurbulentKineticEnergy);
    }

    const auto& r_node = mrData.GetGeometry()[ConditionNodeIndex];
    const double epsilon =
        r_node.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);

    if (epsilon > 0.0) {
        const double tke = r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
        nu_t_derivative = mrData.mCmu * 2.0 * rN[ConditionNodeIndex] * tke / epsilon;
    }

    const double coeff =
        (mrData.mKappa * std::pow(mrData.mYPlus * mrData.mKinematicViscosity, 2));

    return (mrData.mKinematicViscosity +
            mrData.mTurbulentKinematicViscosity / mrData.mEpsilonSigma) *
               5.0 * std::pow(mrData.mUTau, 4) * u_tau_derivative / coeff +
           std::pow(mrData.mUTau, 5) * nu_t_derivative / (coeff * mrData.mEpsilonSigma);
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
    // here we compute derivatives w.r.t. condition nodes
    const auto& r_node = mrData.GetGeometry()[ConditionNodeIndex];
    const double epsilon = r_node.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
    double nu_t_derivative = 0.0;
    if (epsilon > 0.0) {
        const double tke = r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
        nu_t_derivative =
            mrData.mCmu * rN[ConditionNodeIndex] * (-1.0) * std::pow(tke / epsilon, 2);
    }

    const double coeff = (mrData.mKappa * std::pow(mrData.mYPlus * mrData.mKinematicViscosity, 2));

    return std::pow(mrData.mUTau, 5) * nu_t_derivative / (coeff * mrData.mEpsilonSigma);
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

        double u_tau_derivative;
        RansAdjointUtilities::CalculateYPlusAndUtauDerivative(
            y_plus_derivative, u_tau_derivative, mrData.mYPlus,
            mrData.mWallVelocityMagnitude, 0.0, mrData.mWallHeight, wall_height_derivative,
            mrData.mKinematicViscosity, mrData.mKappa, mrData.mBeta, mrData.mYPlusLimit);

        // this does not have any effect since mNumberOfGaussPoints is 1.
        // otherwise this may result in undesired behaviour
    }

    return (mrData.mKinematicViscosity +
            mrData.mTurbulentKinematicViscosity / mrData.mEpsilonSigma) *
           std::pow(mrData.mUTau, 5) *
           (-2.0 * y_plus_derivative / std::pow(mrData.mYPlus, 3)) /
           (mrData.mKappa * std::pow(mrData.mKinematicViscosity, 2));
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
        const double wall_height_derivative =
            RansAdjointUtilities::CalculateWallHeightParentElementDerivative(
                mrData.GetGeometry(), mrData.mrParentElementGeometry,
                DirectionIndex, mrData.mUnitNormal, ZeroVector(3));

        double u_tau_derivative;
        RansAdjointUtilities::CalculateYPlusAndUtauDerivative(
            y_plus_derivative, u_tau_derivative, mrData.mYPlus,
            mrData.mWallVelocityMagnitude, 0.0, mrData.mWallHeight, wall_height_derivative,
            mrData.mKinematicViscosity, mrData.mKappa, mrData.mBeta, mrData.mYPlusLimit);

        // this does not have any effect since mNumberOfGaussPoints is 1.
        // otherwise this may result in undesired behaviour
    }

    return (mrData.mKinematicViscosity +
            mrData.mTurbulentKinematicViscosity / mrData.mEpsilonSigma) *
           std::pow(mrData.mUTau, 5) *
           (-2.0 * y_plus_derivative / std::pow(mrData.mYPlus, 3)) /
           (mrData.mKappa * std::pow(mrData.mKinematicViscosity, 2));
}

// template instantiations

template class EpsilonKBasedWallConditionDataDerivatives<2>;
template class EpsilonKBasedWallConditionDataDerivatives<3>;

} // namespace KEpsilonWallConditionData

} // namespace Kratos