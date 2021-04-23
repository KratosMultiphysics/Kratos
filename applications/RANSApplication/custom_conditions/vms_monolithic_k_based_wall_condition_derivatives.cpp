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
#include "containers/variable.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"

// Application includes
#include "custom_utilities/rans_adjoint_utilities.h"
#include "custom_utilities/fluid_calculation_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Derivative data type includes
#include "custom_conditions/vms_monolithic_k_based_wall_condition_derivative_utilities.h"

// Include base h
#include "vms_monolithic_k_based_wall_condition_derivatives.h"

namespace Kratos
{

/***************************************************************************************************/
/***************************************** Static Methods ******************************************/
/***************************************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes>
void VMSMonolithicKBasedWallConditionDerivatives<TDim, TNumNodes>::Check(
    const Condition& rCondition,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rProcessInfo.Has(TURBULENCE_RANS_C_MU))
        << "TURBULENCE_RANS_C_MU is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rProcessInfo.Has(VON_KARMAN))
        << "VON_KARMAN is not found in process info.\n";

    KRATOS_ERROR_IF_NOT(rCondition.Has(NEIGHBOUR_ELEMENTS))
        << "NEIGHBOUR_ELEMENTS is not found for condition. [ Condition.Id() = "
        << rCondition.Id() << " ].\n";
    KRATOS_ERROR_IF(rCondition.GetValue(NEIGHBOUR_ELEMENTS).size() != 1)
        << "There can be only one parent element for the condition in "
        "NEIGHBOUR_ELEMENTS variable. [ Condition.Id() = "
        << rCondition.Id() << ", NEIGHBOUR_ELEMENTS.size() = "
        << rCondition.GetValue(NEIGHBOUR_ELEMENTS).size() << " ].\n";

    const auto& r_element_properties =
        rCondition.GetValue(NEIGHBOUR_ELEMENTS)[0].GetProperties();
    KRATOS_ERROR_IF_NOT(r_element_properties.Has(DENSITY))
        << "DENSITY is not found in parent element properties. [ "
        "Condition.Id() = "
        << rCondition.Id()
        << ", Properties.Id() = " << r_element_properties.Id() << " ].\n";
    KRATOS_ERROR_IF_NOT(r_element_properties.Has(DYNAMIC_VISCOSITY))
        << "DYNAMIC_VISCOSITY is not found in parent element properties. [ "
        "Condition.Id() = "
        << rCondition.Id()
        << ", Properties.Id() = " << r_element_properties.Id() << " ].\n";

    const auto& r_condition_properties = rCondition.GetProperties();
    KRATOS_ERROR_IF_NOT(r_condition_properties.Has(WALL_SMOOTHNESS_BETA))
        << "WALL_SMOOTHNESS_BETA is not found in condition properties. [ "
        "Condition.Id() = "
        << rCondition.Id()
        << ", Properties.Id() = " << r_condition_properties.Id() << " ].\n";
    KRATOS_ERROR_IF_NOT(r_condition_properties.Has(RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT))
        << "RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT is not found in condition properties. [ "
        "Condition.Id() = "
        << rCondition.Id()
        << ", Properties.Id() = " << r_condition_properties.Id() << " ].\n";

    KRATOS_ERROR_IF_NOT(rCondition.Has(NORMAL))
        << "NORMAL is not found in condition. [ Condition.Id() = " << rCondition.Id()
        << " ].\n";

    KRATOS_ERROR_IF(norm_2(rCondition.GetValue(NORMAL)) == 0.0)
        << "NORMAL is not properly defined for condition. [ Condition.Id() = "
        << rCondition.Id() << " ].\n";

    KRATOS_ERROR_IF_NOT(rCondition.Has(NORMAL_SHAPE_DERIVATIVE))
        << "NORMAL_SHAPE_DERIVATIVE is not found in condition. [ "
        "Condition.Id() = "
        << rCondition.Id() << " ].\n";

    KRATOS_ERROR_IF_NOT(rCondition.Has(DISTANCE))
        << "DISTANCE is not found in condition. [ "
        "Condition.Id() = "
        << rCondition.Id() << " ].\n";

    const auto& r_geometry = rCondition.GetGeometry();
    for (IndexType i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node) {
        const auto& r_node = r_geometry[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
GeometryData::IntegrationMethod VMSMonolithicKBasedWallConditionDerivatives<TDim, TNumNodes>::GetIntegrationMethod()
{
    return GeometryData::IntegrationMethod::GI_GAUSS_2;
}

template <unsigned int TDim, unsigned int TNumNodes>
void VMSMonolithicKBasedWallConditionDerivatives<TDim, TNumNodes>::InitializeCondition(
    Condition& rCondition,
    const ProcessInfo& rProcessInfo)
{
    const double wall_height = RansCalculationUtilities::CalculateWallHeight(rCondition, rCondition.GetValue(NORMAL));
    rCondition.SetValue(DISTANCE, wall_height);

    const IndexType number_of_gauss_points = rCondition.GetGeometry().IntegrationPointsNumber(GeometryData::IntegrationMethod::GI_GAUSS_2);
    if (!rCondition.Has(GAUSS_RANS_Y_PLUS) || rCondition.GetValue(GAUSS_RANS_Y_PLUS).size() != number_of_gauss_points) {
        rCondition.SetValue(GAUSS_RANS_Y_PLUS, Vector(number_of_gauss_points));
    }
}

/***************************************************************************************************/
/********************************************** Data ***********************************************/
/***************************************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes>
VMSMonolithicKBasedWallConditionDerivatives<TDim, TNumNodes>::Data::Data(
    Condition& rCondition,
    Element& rParentElement,
    const ProcessInfo& rProcessInfo)
    : mrCondition(rCondition),
      mrParentElement(rParentElement)
{
    KRATOS_TRY

    mGaussIntegrationPointIndex = 0;

    mCmu = rProcessInfo[TURBULENCE_RANS_C_MU];
    mKappa = rProcessInfo[VON_KARMAN];

    const auto& element_properties = rParentElement.GetProperties();
    mDensity = element_properties[DENSITY];
    mKinematicViscosity = element_properties[DYNAMIC_VISCOSITY] / mDensity;

    const auto& condition_properties = rCondition.GetProperties();
    mBeta = condition_properties[WALL_SMOOTHNESS_BETA];
    mYPlusLimit = condition_properties[RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT];

    const auto& normal = rCondition.GetValue(NORMAL);
    noalias(mUnitNormal) = normal / norm_2(normal);
    mNormalMagnitude = rCondition.GetGeometry().Area();

    mCmu25 = std::pow(mCmu, 0.25);
    mInvKappa = 1.0 / mKappa;

    mWallHeight = rCondition.GetValue(DISTANCE);

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void VMSMonolithicKBasedWallConditionDerivatives<TDim, TNumNodes>::Data::CalculateGaussPointData(
    const double W,
    const Vector& rN,
    const int Step)
{
    KRATOS_TRY

    FluidCalculationUtilities::EvaluateInPoint(
        mrCondition.GetGeometry(), rN, Step,
        std::tie(mTurbulentKineticEnergy, TURBULENT_KINETIC_ENERGY),
        std::tie(mWallVelocity, VELOCITY));

    mWallVelocityMagnitude = norm_2(mWallVelocity);

    RansCalculationUtilities::CalculateYPlusAndUtau(
        mYPlus, mUTau, mWallVelocityMagnitude, mWallHeight, mKinematicViscosity, mKappa, mBeta);
    mYPlus = std::max(mYPlus, mYPlusLimit);

    double& y_plus = mrCondition.GetValue(GAUSS_RANS_Y_PLUS)[mGaussIntegrationPointIndex++];
    y_plus = mYPlus;

    mTkeBasedUTau = mCmu25 * std::sqrt(std::max(mTurbulentKineticEnergy, 0.0));
    mUBasedUTau = mWallVelocityMagnitude / (mInvKappa * std::log(mYPlus) + mBeta);

    mUTau = std::max(mTkeBasedUTau, mUBasedUTau);

    KRATOS_CATCH("");
}

/***************************************************************************************************/
/************************************* Residual Contributions **************************************/
/***************************************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes>
VMSMonolithicKBasedWallConditionDerivatives<TDim, TNumNodes>::ResidualsContributions::ResidualsContributions(
    Data& rData)
    : mrData(rData)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
void VMSMonolithicKBasedWallConditionDerivatives<TDim, TNumNodes>::ResidualsContributions::AddGaussPointResidualsContributions(
    VectorF& rResidual,
    const double W,
    const Vector& rN)
{
    KRATOS_TRY

    if (mrData.mWallVelocityMagnitude > 1e-16) {
        const double value = mrData.mDensity * std::pow(mrData.mUTau, 2) * W /
                             mrData.mWallVelocityMagnitude;

        for (IndexType a = 0; a < TNumNodes; ++a) {
            for (IndexType i = 0; i < TDim; ++i) {
                rResidual[a * TBlockSize + i] -= rN[a] * value * mrData.mWallVelocity[i];
            }
        }
    }

    KRATOS_CATCH("");
}

/***************************************************************************************************/
/************************************** Variable Derivatives ***************************************/
/***************************************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes>
template <class TDerivativesType>
VMSMonolithicKBasedWallConditionDerivatives<TDim, TNumNodes>::VariableDerivatives<
    TDerivativesType>::VariableDerivatives(Data& rData)
    : mrData(rData),
      mResidualWeightDerivativeContributions(rData)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
template <class TDerivativesType>
void VMSMonolithicKBasedWallConditionDerivatives<TDim, TNumNodes>::VariableDerivatives<TDerivativesType>::CalculateGaussPointResidualsDerivativeContributions(
    VectorF& rResidualDerivative,
    const IndexType ConditionNodeIndex,
    const IndexType ParentElementNodeIndex,
    const IndexType DirectionIndex,
    const double W,
    const Vector& rN,
    const double WDerivative,
    const double DetJDerivative)
{
    KRATOS_TRY

    rResidualDerivative.clear();

    if (mrData.mWallVelocityMagnitude > 1e-16) {
        const auto& wall_velocity_derivative = TDerivativesType::CalculateWallVelocityDerivative(
            ConditionNodeIndex, DirectionIndex, rN);

        const auto& wall_height_derivative = TDerivativesType::CalculateWallHeightConditionDerivative(
            mrData.mrCondition.GetGeometry(), mrData.mrParentElement.GetGeometry(),
            ConditionNodeIndex, DirectionIndex, mrData.mNormalMagnitude, mrData.mUnitNormal);

        const double wall_velocity_magnitude_derivative =
            RansAdjointUtilities::CalculateVectorNormDerivative(
                mrData.mWallVelocityMagnitude, mrData.mWallVelocity, wall_velocity_derivative);

        double y_plus_derivative, u_tau_derivative;
        RansAdjointUtilities::CalculateYPlusAndUtauDerivative(
            y_plus_derivative, u_tau_derivative, mrData.mYPlus,
            mrData.mWallVelocityMagnitude, wall_velocity_magnitude_derivative,
            mrData.mWallHeight, wall_height_derivative, mrData.mKinematicViscosity,
            mrData.mKappa, mrData.mBeta, mrData.mYPlusLimit);

        y_plus_derivative = (mrData.mYPlus > mrData.mYPlusLimit) ? y_plus_derivative : 0.0;

        if (mrData.mTkeBasedUTau > mrData.mUBasedUTau) {
            if (mrData.mTurbulentKineticEnergy > 0.0) {
                u_tau_derivative = mrData.mCmu25 * 0.5 * rN[ConditionNodeIndex] * TurbulentKineticEnergyDerivativeFactor / std::sqrt(mrData.mTurbulentKineticEnergy);
            } else {
                u_tau_derivative = 0.0;
            }
        } else {
            const double coeff1 = (mrData.mInvKappa * std::log(mrData.mYPlus) + mrData.mBeta);
            u_tau_derivative = wall_velocity_magnitude_derivative / coeff1 - mrData.mWallVelocityMagnitude * (mrData.mInvKappa * y_plus_derivative / mrData.mYPlus) / std::pow(coeff1, 2);
        }

        const double value = mrData.mDensity * std::pow(mrData.mUTau, 2) * W / mrData.mWallVelocityMagnitude;
        const double value_derivative = mrData.mDensity * W * (2.0 * mrData.mUTau * u_tau_derivative / mrData.mWallVelocityMagnitude - std::pow(mrData.mUTau / mrData.mWallVelocityMagnitude, 2) * wall_velocity_magnitude_derivative);

        for (IndexType a = 0; a < TNumNodes; ++a) {
            rResidualDerivative[a * TBlockSize + DirectionIndex] -= rN[a] * value * rN[ConditionNodeIndex] * VelocityDerivativeFactor;
            for (IndexType i = 0; i < TDim; ++i) {
                rResidualDerivative[a * TBlockSize + i] -= rN[a] * value_derivative * mrData.mWallVelocity[i];
            }
        }

        mResidualWeightDerivativeContributions.AddGaussPointResidualsContributions(rResidualDerivative, WDerivative, rN);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
template <class TDerivativesType>
void VMSMonolithicKBasedWallConditionDerivatives<TDim, TNumNodes>::VariableDerivatives<TDerivativesType>::CalculateGaussPointResidualsDerivativeContributions(
    VectorF& rResidualDerivative,
    const IndexType ParentElementNodeIndex,
    const IndexType DirectionIndex,
    const double W,
    const Vector& rN)
{
    KRATOS_TRY

    rResidualDerivative.clear();

    if (mrData.mWallVelocityMagnitude > 1e-16) {
        const auto& wall_height_derivative = TDerivativesType::CalculateWallHeightParentElementDerivative(
            mrData.mrCondition.GetGeometry(), mrData.mrParentElement.GetGeometry(),
            DirectionIndex, mrData.mUnitNormal);

        double y_plus_derivative, u_tau_derivative;
        RansAdjointUtilities::CalculateYPlusAndUtauDerivative(
            y_plus_derivative, u_tau_derivative, mrData.mYPlus,
            mrData.mWallVelocityMagnitude, 0.0,
            mrData.mWallHeight, wall_height_derivative, mrData.mKinematicViscosity,
            mrData.mKappa, mrData.mBeta, mrData.mYPlusLimit);

        y_plus_derivative = (mrData.mYPlus > mrData.mYPlusLimit) ? y_plus_derivative : 0.0;

        if (mrData.mTkeBasedUTau > mrData.mUBasedUTau) {
            u_tau_derivative = 0.0;
        } else {
            u_tau_derivative = -mrData.mWallVelocityMagnitude * (mrData.mInvKappa * y_plus_derivative / mrData.mYPlus) / std::pow(mrData.mInvKappa * std::log(mrData.mYPlus) + mrData.mBeta, 2);
        }

        const double value_derivative = mrData.mDensity * W * 2.0 * mrData.mUTau * u_tau_derivative / mrData.mWallVelocityMagnitude;

        for (IndexType a = 0; a < TNumNodes; ++a) {
            for (IndexType i = 0; i < TDim; ++i) {
                rResidualDerivative[a * TBlockSize + i] -= rN[a] * value_derivative * mrData.mWallVelocity[i];
            }
        }
    }

    KRATOS_CATCH("");
}

// template instantiations

template class VMSMonolithicKBasedWallConditionDerivatives<2, 2>;
template class VMSMonolithicKBasedWallConditionDerivatives<2, 2>::VariableDerivatives<typename VMSMonolithicKBasedWallConditionDerivativeUtilities<2>::VelocityDerivative>;
template class VMSMonolithicKBasedWallConditionDerivatives<2, 2>::VariableDerivatives<typename VMSMonolithicKBasedWallConditionDerivativeUtilities<2>::KDerivative>;
template class VMSMonolithicKBasedWallConditionDerivatives<2, 2>::VariableDerivatives<typename VMSMonolithicKBasedWallConditionDerivativeUtilities<2>::ShapeDerivative>;
template class VMSMonolithicKBasedWallConditionDerivatives<2, 2>::VariableDerivatives<typename VMSMonolithicKBasedWallConditionDerivativeUtilities<2>::NonRelatedDerivative>;


template class VMSMonolithicKBasedWallConditionDerivatives<3, 3>;
template class VMSMonolithicKBasedWallConditionDerivatives<3, 3>::VariableDerivatives<typename VMSMonolithicKBasedWallConditionDerivativeUtilities<3>::VelocityDerivative>;
template class VMSMonolithicKBasedWallConditionDerivatives<3, 3>::VariableDerivatives<typename VMSMonolithicKBasedWallConditionDerivativeUtilities<3>::KDerivative>;
template class VMSMonolithicKBasedWallConditionDerivatives<3, 3>::VariableDerivatives<typename VMSMonolithicKBasedWallConditionDerivativeUtilities<3>::ShapeDerivative>;
template class VMSMonolithicKBasedWallConditionDerivatives<3, 3>::VariableDerivatives<typename VMSMonolithicKBasedWallConditionDerivativeUtilities<3>::NonRelatedDerivative>;

} // namespace Kratos
