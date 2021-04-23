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
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

// Application includes
#include "custom_elements/data_containers/k_epsilon/element_data_derivative_utilities.h"
#include "custom_elements/data_containers/k_epsilon/element_data_utilities.h"
#include "custom_elements/data_containers/k_omega_sst/element_data_derivative_utilities.h"
#include "custom_elements/data_containers/k_omega_sst/element_data_utilities.h"
#include "custom_utilities/fluid_calculation_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "k_element_data_derivatives.h"

namespace Kratos
{

namespace KOmegaSSTElementData
{

/***************************************************************/
/************************ Element Data *************************/
/***************************************************************/

template <unsigned int TDim, unsigned int TNumNodes>
const Variable<double>& KElementDataDerivatives<TDim, TNumNodes>::Data::GetAdjointScalarVariable()
{
    return RANS_SCALAR_1_ADJOINT_1;
}

template <unsigned int TDim, unsigned int TNumNodes>
void KElementDataDerivatives<TDim, TNumNodes>::Data::Check(
    const Element& rElement,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY


    const auto& r_geometry = rElement.GetGeometry();
    const auto& r_properties = rElement.GetProperties();
    const int number_of_nodes = r_geometry.PointsNumber();

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENT_KINETIC_ENERGY_SIGMA_1))
        << "TURBULENT_KINETIC_ENERGY_SIGMA_1 is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENT_KINETIC_ENERGY_SIGMA_2))
        << "TURBULENT_KINETIC_ENERGY_SIGMA_2 is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_2))
        << "TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_2 is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENCE_RANS_C_MU))
        << "TURBULENCE_RANS_C_MU is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENCE_RANS_A1))
        << "TURBULENCE_RANS_A1 is not found in process info.\n";

    KRATOS_ERROR_IF_NOT(r_properties.Has(DYNAMIC_VISCOSITY))
        << "DYNAMIC_VISCOSITY is not found in element properties [ Element.Id() = "
        << rElement.Id() << ", Properties.Id() = " << r_properties.Id() << " ].\n";
    KRATOS_ERROR_IF_NOT(r_properties.Has(DENSITY))
        << "DENSITY is not found in element properties [ Element.Id() = "
        << rElement.Id() << ", Properties.Id() = " << r_properties.Id() << " ].\n";

    for (int i_node = 0; i_node < number_of_nodes; ++i_node) {
        const auto& r_node = r_geometry[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_AUXILIARY_VARIABLE_1, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_SCALAR_1_ADJOINT_1, r_node);

        KRATOS_CHECK_DOF_IN_NODE(RANS_SCALAR_1_ADJOINT_1, r_node);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void KElementDataDerivatives<TDim, TNumNodes>::Data::CalculateGaussPointData(
    const Vector& rN,
    const Matrix& rdNdX,
    const int Step)
{
    using namespace KOmegaSSTElementData;

    const auto& r_geometry = this->GetGeometry();

    FluidCalculationUtilities::EvaluateInPoint(
        r_geometry, rN, Step,
        std::tie(mTurbulentKineticEnergy, TURBULENT_KINETIC_ENERGY),
        std::tie(mTurbulentSpecificEnergyDissipationRate, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE),
        std::tie(mWallDistance, DISTANCE),
        std::tie(mEffectiveVelocity, VELOCITY));

    FluidCalculationUtilities::EvaluateGradientInPoint(
        this->GetGeometry(), rdNdX,
        std::tie(mTurbulentKineticEnergyGradient, TURBULENT_KINETIC_ENERGY),
        std::tie(mTurbulentSpecificEnergyDissipationRateGradient, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE),
        std::tie(mVelocityGradient, VELOCITY));

    mCrossDiffusion = KOmegaSSTElementData::CalculateCrossDiffusionTerm<TDim>(
        mSigmaOmega2, mTurbulentSpecificEnergyDissipationRate,
        mTurbulentKineticEnergyGradient, mTurbulentSpecificEnergyDissipationRateGradient);

    mF1 = KOmegaSSTElementData::CalculateF1(
        mTurbulentKineticEnergy, mTurbulentSpecificEnergyDissipationRate,
        mKinematicViscosity, mWallDistance, mBetaStar, mCrossDiffusion, mSigmaOmega2);

    mBlendedSimgaK = KOmegaSSTElementData::CalculateBlendedPhi(mSigmaK1, mSigmaK2, mF1);

    mVelocityDivergence = RansCalculationUtilities::CalculateMatrixTrace<TDim>(mVelocityGradient);

    mF2 = KOmegaSSTElementData::CalculateF2(mTurbulentKineticEnergy, mTurbulentSpecificEnergyDissipationRate, mKinematicViscosity, mWallDistance, mBetaStar);

    mSymmetricVelocityGradient = (mVelocityGradient + trans(mVelocityGradient)) * 0.5;
    mT = norm_frobenius(mSymmetricVelocityGradient) * 1.414;

    mTurbulentKinematicViscosity = KOmegaSSTElementData::CalculateTurbulentKinematicViscosity(mTurbulentKineticEnergy, mTurbulentSpecificEnergyDissipationRate, mT, mF2, mA1);
    mTurbulentKinematicViscosity = std::max(mTurbulentKinematicViscosity, 1e-12);

    for (IndexType i = 0; i < TNumNodes; ++i) {
        const auto& velocity = r_geometry[i].FastGetSolutionStepValue(VELOCITY, Step);
        for (IndexType j = 0; j < TDim; ++j) {
            mNodalVelocity(i, j) = velocity[j];
        }
    }

    mEffectiveKinematicViscosity = mKinematicViscosity + mBlendedSimgaK * mTurbulentKinematicViscosity;
    mReactionTerm = std::max(mBetaStar * mTurbulentKineticEnergy / mTurbulentKinematicViscosity + (2.0 / 3.0) * mVelocityDivergence, 0.0);
    mSourceTerm = KEpsilonElementData::CalculateProductionTerm<TDim>(mVelocityGradient, mTurbulentKinematicViscosity);
}

/***************************************************************/
/********************* Velocity Derivative *********************/
/***************************************************************/

template <unsigned int TDim, unsigned int TNumNodes>
KElementDataDerivatives<TDim, TNumNodes>::UDerivative::UDerivative(
    const DataType& rData,
    const IndexType NodeIndex,
    const IndexType DirectionIndex,
    const double W,
    const Vector& rN,
    const Matrix& rdNdX,
    const double WDerivative,
    const double DetJDerivative,
    const Matrix& rdNdXDerivative)
    : BaseType(
        NodeIndex,
        DirectionIndex,
        rData.GetGeometry(),
        W,
        rN,
        rdNdX,
        WDerivative,
        DetJDerivative,
        rdNdXDerivative),
        mrData(rData)
{
    MatrixDD symmetric_velocity_gradient_derivative = ZeroMatrix(TDim, TDim);

    for (IndexType i = 0; i < TDim; ++i) {
        symmetric_velocity_gradient_derivative(this->mDirectionIndex, i) += 0.5 * this->mrdNdX(this->mNodeIndex, i);
        symmetric_velocity_gradient_derivative(i, this->mDirectionIndex) += 0.5 * this->mrdNdX(this->mNodeIndex, i);
    }

    double t_derivative = 0.0;
    for (IndexType i = 0; i < TDim; ++i) {
        for (IndexType j = 0; j < TDim; ++j) {
            t_derivative += this->mrData.mSymmetricVelocityGradient(i, j) * symmetric_velocity_gradient_derivative(i, j);
        }
    }

    t_derivative *= (1.414 * 1.414 / this->mrData.mT);

    const double f2_derivative = KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateF2Derivative(
        mrData.mTurbulentKineticEnergy, 0.0, mrData.mTurbulentSpecificEnergyDissipationRate,
        0.0, mrData.mKinematicViscosity, mrData.mWallDistance, 0.0, mrData.mBetaStar);

    mGaussTurbulentKinematicViscosityDerivative = 0.0;
    if (mrData.mTurbulentKinematicViscosity > 1e-12) {
        mGaussTurbulentKinematicViscosityDerivative =
            KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateTurbulentKinematicViscosityDerivative(
                mrData.mTurbulentKineticEnergy, 0.0,
                mrData.mTurbulentSpecificEnergyDissipationRate, 0.0, mrData.mT,
                t_derivative, mrData.mF2, f2_derivative, mrData.mA1);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
const Variable<double>& KElementDataDerivatives<TDim, TNumNodes>::UDerivative::GetDerivativeVariable() const
{
    switch (this->mDirectionIndex) {
    case 0:
        return VELOCITY_X;
        break;
    case 1:
        return VELOCITY_Y;
        break;
    case 2:
        return VELOCITY_Z;
        break;
    default:
        return Variable<double>::StaticObject();
    };
}

template <unsigned int TDim, unsigned int TNumNodes>
array_1d<double, TDim> KElementDataDerivatives<TDim, TNumNodes>::UDerivative::CalculateEffectiveVelocityDerivative() const
{
    array_1d<double, TDim> result = ZeroVector(TDim);
    result[this->mDirectionIndex] = this->mrN[this->mNodeIndex];
    return result;
}

template <unsigned int TDim, unsigned int TNumNodes>
double KElementDataDerivatives<TDim, TNumNodes>::UDerivative::CalculateEffectiveKinematicViscosityDerivative() const
{
    const ArrayD& zero = ZeroVector(TDim);

    const double cross_diffusion_derivative =
        KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateCrossDiffusionTermDerivative(
            mrData.mSigmaOmega2, mrData.mTurbulentSpecificEnergyDissipationRate,
            0.0, mrData.mTurbulentKineticEnergyGradient, zero,
            mrData.mTurbulentSpecificEnergyDissipationRateGradient, zero);

    const double f1_derivative = KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateF1Derivative(
        mrData.mTurbulentKineticEnergy, 0.0, mrData.mTurbulentSpecificEnergyDissipationRate,
        0.0, mrData.mKinematicViscosity, mrData.mWallDistance, 0.0, mrData.mBetaStar,
        mrData.mCrossDiffusion, cross_diffusion_derivative, mrData.mSigmaOmega2);

    const double blended_sigma_k_derivative =
        KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateBlendedPhiDerivative(
            mrData.mSigmaK1, 0.0, mrData.mSigmaK2, 0.0, mrData.mF1, f1_derivative);

    return blended_sigma_k_derivative * mrData.mTurbulentKinematicViscosity +
           mrData.mBlendedSimgaK * mGaussTurbulentKinematicViscosityDerivative;
}

template <unsigned int TDim, unsigned int TNumNodes>
double KElementDataDerivatives<TDim, TNumNodes>::UDerivative::CalculateReactionTermDerivative() const
{
    if (mrData.mReactionTerm > 0.0) {
        return this->mrdNdX(this->mNodeIndex, this->mDirectionIndex) * (2.0 / 3.0) -
               mrData.mBetaStar * mrData.mTurbulentKineticEnergy *
                   mGaussTurbulentKinematicViscosityDerivative /
                   std::pow(mrData.mTurbulentKinematicViscosity, 2);
    } else {
        return 0.0;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
double KElementDataDerivatives<TDim, TNumNodes>::UDerivative::CalculateSourceTermDerivative() const
{
    return KEpsilonElementData::AdjointUtilities<TDim, TNumNodes>::CalculateProductionVelocityDerivative(
        this->mNodeIndex, this->mDirectionIndex, this->mrData.mSourceTerm,
        this->mrData.mTurbulentKinematicViscosity,
        mGaussTurbulentKinematicViscosityDerivative,
        this->mrData.mVelocityGradient, this->mrdNdX);
}

/***************************************************************/
/************************ K Derivative  ************************/
/***************************************************************/

template <unsigned int TDim, unsigned int TNumNodes>
KElementDataDerivatives<TDim, TNumNodes>::KDerivative::KDerivative(
    const Data& rData,
    const IndexType NodeIndex,
    const IndexType DirectionIndex,
    const double W,
    const Vector& rN,
    const Matrix& rdNdX,
    const double WDerivative,
    const double DetJDerivative,
    const Matrix& rdNdXDerivative)
    : BaseType(NodeIndex, DirectionIndex, rData.GetGeometry(), W, rN, rdNdX, WDerivative, DetJDerivative, rdNdXDerivative),
      mrData(rData)
{
    const double tke_derivative = this->mrN[this->mNodeIndex];

    const double f2_derivative = KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateF2Derivative(
        mrData.mTurbulentKineticEnergy, tke_derivative, mrData.mTurbulentSpecificEnergyDissipationRate,
        0.0, mrData.mKinematicViscosity, mrData.mWallDistance, 0.0, mrData.mBetaStar);

    mGaussTurbulentKinematicViscosityDerivative = 0.0;
    if (mrData.mTurbulentKinematicViscosity > 1e-12) {
        mGaussTurbulentKinematicViscosityDerivative =
            KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateTurbulentKinematicViscosityDerivative(
                mrData.mTurbulentKineticEnergy, tke_derivative,
                mrData.mTurbulentSpecificEnergyDissipationRate, 0.0, mrData.mT,
                0.0, mrData.mF2, f2_derivative, mrData.mA1);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
const Variable<double>& KElementDataDerivatives<TDim, TNumNodes>::KDerivative::GetDerivativeVariable() const
{
    return TURBULENT_KINETIC_ENERGY;
}

template <unsigned int TDim, unsigned int TNumNodes>
array_1d<double, TDim> KElementDataDerivatives<TDim, TNumNodes>::KDerivative::CalculateEffectiveVelocityDerivative() const
{
    return ZeroVector(TDim);
}

template <unsigned int TDim, unsigned int TNumNodes>
double KElementDataDerivatives<TDim, TNumNodes>::KDerivative::CalculateEffectiveKinematicViscosityDerivative() const
{
    const ArrayD& zero = ZeroVector(TDim);

    const double cross_diffusion_derivative =
        KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateCrossDiffusionTermDerivative(
            mrData.mSigmaOmega2, mrData.mTurbulentSpecificEnergyDissipationRate,
            0.0, mrData.mTurbulentKineticEnergyGradient, row(this->mrdNdX, this->mNodeIndex),
            mrData.mTurbulentSpecificEnergyDissipationRateGradient, zero);

    const double f1_derivative = KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateF1Derivative(
        mrData.mTurbulentKineticEnergy, this->mrN[this->mNodeIndex], mrData.mTurbulentSpecificEnergyDissipationRate,
        0.0, mrData.mKinematicViscosity, mrData.mWallDistance, 0.0, mrData.mBetaStar,
        mrData.mCrossDiffusion, cross_diffusion_derivative, mrData.mSigmaOmega2);

    const double blended_sigma_k_derivative =
        KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateBlendedPhiDerivative(
            mrData.mSigmaK1, 0.0, mrData.mSigmaK2, 0.0, mrData.mF1, f1_derivative);

    return blended_sigma_k_derivative * mrData.mTurbulentKinematicViscosity +
           mrData.mBlendedSimgaK * mGaussTurbulentKinematicViscosityDerivative;
}

template <unsigned int TDim, unsigned int TNumNodes>
double KElementDataDerivatives<TDim, TNumNodes>::KDerivative::CalculateReactionTermDerivative() const
{
    if (mrData.mReactionTerm > 0.0) {
        double value = 0.0;
        value += this->mrN[this->mNodeIndex] / mrData.mTurbulentKinematicViscosity;
        value -= mrData.mTurbulentKineticEnergy *  mGaussTurbulentKinematicViscosityDerivative / std::pow(mrData.mTurbulentKinematicViscosity, 2);

        return value * mrData.mBetaStar;
    } else {
        return 0.0;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
double KElementDataDerivatives<TDim, TNumNodes>::KDerivative::CalculateSourceTermDerivative() const
{
    return KEpsilonElementData::AdjointUtilities<TDim, TNumNodes>::CalculateProductionScalarDerivative(
        mrData.mTurbulentKinematicViscosity, mrData.mSourceTerm,
        mGaussTurbulentKinematicViscosityDerivative);
}

/***************************************************************/
/********************** Omega Derivative *********************/
/***************************************************************/

template <unsigned int TDim, unsigned int TNumNodes>
KElementDataDerivatives<TDim, TNumNodes>::OmegaDerivative::OmegaDerivative(
    const Data& rData,
    const IndexType NodeIndex,
    const IndexType DirectionIndex,
    const double W,
    const Vector& rN,
    const Matrix& rdNdX,
    const double WDerivative,
    const double DetJDerivative,
    const Matrix& rdNdXDerivative)
    : BaseType(NodeIndex, DirectionIndex, rData.GetGeometry(), W, rN, rdNdX, WDerivative, DetJDerivative, rdNdXDerivative),
      mrData(rData)
{
    const double omega_derivative = this->mrN[this->mNodeIndex];

    const double f2_derivative = KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateF2Derivative(
        mrData.mTurbulentKineticEnergy, 0.0, mrData.mTurbulentSpecificEnergyDissipationRate,
        omega_derivative, mrData.mKinematicViscosity, mrData.mWallDistance, 0.0, mrData.mBetaStar);

    mGaussTurbulentKinematicViscosityDerivative = 0.0;
    if (mrData.mTurbulentKinematicViscosity > 1e-12) {
        mGaussTurbulentKinematicViscosityDerivative =
            KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateTurbulentKinematicViscosityDerivative(
                mrData.mTurbulentKineticEnergy, 0.0,
                mrData.mTurbulentSpecificEnergyDissipationRate, omega_derivative, mrData.mT,
                0.0, mrData.mF2, f2_derivative, mrData.mA1);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
const Variable<double>& KElementDataDerivatives<TDim, TNumNodes>::OmegaDerivative::GetDerivativeVariable() const
{
    return TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE;
}

template <unsigned int TDim, unsigned int TNumNodes>
array_1d<double, TDim> KElementDataDerivatives<TDim, TNumNodes>::OmegaDerivative::CalculateEffectiveVelocityDerivative() const
{
    return ZeroVector(TDim);
}

template <unsigned int TDim, unsigned int TNumNodes>
double KElementDataDerivatives<TDim, TNumNodes>::OmegaDerivative::CalculateEffectiveKinematicViscosityDerivative() const
{
    const ArrayD& zero = ZeroVector(TDim);

    const double cross_diffusion_derivative =
        KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateCrossDiffusionTermDerivative(
            mrData.mSigmaOmega2, mrData.mTurbulentSpecificEnergyDissipationRate,
            this->mrN[this->mNodeIndex], mrData.mTurbulentKineticEnergyGradient, zero,
            mrData.mTurbulentSpecificEnergyDissipationRateGradient, row(this->mrdNdX, this->mNodeIndex));

    const double f1_derivative = KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateF1Derivative(
        mrData.mTurbulentKineticEnergy, 0.0, mrData.mTurbulentSpecificEnergyDissipationRate,
        this->mrN[this->mNodeIndex], mrData.mKinematicViscosity, mrData.mWallDistance, 0.0, mrData.mBetaStar,
        mrData.mCrossDiffusion, cross_diffusion_derivative, mrData.mSigmaOmega2);

    const double blended_sigma_k_derivative =
        KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateBlendedPhiDerivative(
            mrData.mSigmaK1, 0.0, mrData.mSigmaK2, 0.0, mrData.mF1, f1_derivative);

    return blended_sigma_k_derivative * mrData.mTurbulentKinematicViscosity +
           mrData.mBlendedSimgaK * mGaussTurbulentKinematicViscosityDerivative;
}

template <unsigned int TDim, unsigned int TNumNodes>
double KElementDataDerivatives<TDim, TNumNodes>::OmegaDerivative::CalculateReactionTermDerivative() const
{
    if (mrData.mReactionTerm > 0.0) {
        double value = 0.0;
        value -= mrData.mTurbulentKineticEnergy *  mGaussTurbulentKinematicViscosityDerivative / std::pow(mrData.mTurbulentKinematicViscosity, 2);

        return value * mrData.mBetaStar;
    } else {
        return 0.0;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
double KElementDataDerivatives<TDim, TNumNodes>::OmegaDerivative::CalculateSourceTermDerivative() const
{
    return KEpsilonElementData::AdjointUtilities<TDim, TNumNodes>::CalculateProductionScalarDerivative(
        mrData.mTurbulentKinematicViscosity, mrData.mSourceTerm,
        mGaussTurbulentKinematicViscosityDerivative);
}

/***************************************************************/
/*********************** Shape Derivative **********************/
/***************************************************************/

template <unsigned int TDim, unsigned int TNumNodes>
KElementDataDerivatives<TDim, TNumNodes>::ShapeDerivative::ShapeDerivative(
    const DataType& rData,
    const IndexType NodeIndex,
    const IndexType DirectionIndex,
    const double W,
    const Vector& rN,
    const Matrix& rdNdX,
    const double WDerivative,
    const double DetJDerivative,
    const Matrix& rdNdXDerivative)
    : BaseType(
        NodeIndex,
        DirectionIndex,
        rData.GetGeometry(),
        W,
        rN,
        rdNdX,
        WDerivative,
        DetJDerivative,
        rdNdXDerivative),
        mrData(rData)
{
    FluidCalculationUtilities::EvaluateGradientInPoint(
        this->mrData.GetGeometry(), this->mrdNdXDerivative,
        std::tie(mVelocityGradientDerivative, VELOCITY),
        std::tie(mTurbulentKineticEnergyGradientDerivative, TURBULENT_KINETIC_ENERGY),
        std::tie(mTurbulentSpecificEnergyDissipationRateGradientDerivative, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE));

    const MatrixDD& symmetric_velocity_gradient_derivative  = 0.5 * (mVelocityGradientDerivative + trans(mVelocityGradientDerivative));

    double t_derivative = 0.0;
    for (IndexType i = 0; i < TDim; ++i) {
        for (IndexType j = 0; j < TDim; ++j) {
            t_derivative += this->mrData.mSymmetricVelocityGradient(i, j) * symmetric_velocity_gradient_derivative(i, j);
        }
    }

    t_derivative *= (1.414 * 1.414  / this->mrData.mT);

    // TODO: Calculate wall distance derivatives also
    const double f2_derivative = KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateF2Derivative(
        mrData.mTurbulentKineticEnergy, 0.0, mrData.mTurbulentSpecificEnergyDissipationRate,
        0.0, mrData.mKinematicViscosity, mrData.mWallDistance, 0.0, mrData.mBetaStar);

    mGaussTurbulentKinematicViscosityDerivative = 0.0;
    if (mrData.mTurbulentKinematicViscosity > 1e-12) {
        mGaussTurbulentKinematicViscosityDerivative =
            KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateTurbulentKinematicViscosityDerivative(
                mrData.mTurbulentKineticEnergy, 0.0,
                mrData.mTurbulentSpecificEnergyDissipationRate, 0.0, mrData.mT,
                t_derivative, mrData.mF2, f2_derivative, mrData.mA1);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
const Variable<double>& KElementDataDerivatives<TDim, TNumNodes>::ShapeDerivative::GetDerivativeVariable() const
{
    switch (this->mDirectionIndex) {
    case 0:
        return SHAPE_SENSITIVITY_X;
        break;
    case 1:
        return SHAPE_SENSITIVITY_Y;
        break;
    case 2:
        return SHAPE_SENSITIVITY_Z;
        break;
    default:
        return Variable<double>::StaticObject();
    };
}

template <unsigned int TDim, unsigned int TNumNodes>
array_1d<double, TDim> KElementDataDerivatives<TDim, TNumNodes>::ShapeDerivative::CalculateEffectiveVelocityDerivative() const
{
    return ZeroVector(TDim);
}

template <unsigned int TDim, unsigned int TNumNodes>
double KElementDataDerivatives<TDim, TNumNodes>::ShapeDerivative::CalculateEffectiveKinematicViscosityDerivative() const
{
    const double cross_diffusion_derivative =
        KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateCrossDiffusionTermDerivative(
            mrData.mSigmaOmega2, mrData.mTurbulentSpecificEnergyDissipationRate,
            0.0, mrData.mTurbulentKineticEnergyGradient, mTurbulentKineticEnergyGradientDerivative,
            mrData.mTurbulentSpecificEnergyDissipationRateGradient, mTurbulentSpecificEnergyDissipationRateGradientDerivative);

    const double f1_derivative = KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateF1Derivative(
        mrData.mTurbulentKineticEnergy, 0.0, mrData.mTurbulentSpecificEnergyDissipationRate,
        0.0, mrData.mKinematicViscosity, mrData.mWallDistance, 0.0, mrData.mBetaStar,
        mrData.mCrossDiffusion, cross_diffusion_derivative, mrData.mSigmaOmega2);

    const double blended_sigma_k_derivative =
        KOmegaSSTElementData::AdjointUtilities<TDim>::CalculateBlendedPhiDerivative(
            mrData.mSigmaK1, 0.0, mrData.mSigmaK2, 0.0, mrData.mF1, f1_derivative);

    return blended_sigma_k_derivative * mrData.mTurbulentKinematicViscosity +
           mrData.mBlendedSimgaK * mGaussTurbulentKinematicViscosityDerivative;
}

template <unsigned int TDim, unsigned int TNumNodes>
double KElementDataDerivatives<TDim, TNumNodes>::ShapeDerivative::CalculateReactionTermDerivative() const
{
    if (mrData.mReactionTerm > 0.0) {
        return (2.0 / 3.0) * RansCalculationUtilities::CalculateMatrixTrace<TDim>(
                                 mVelocityGradientDerivative) -
               mrData.mBetaStar * mrData.mTurbulentKineticEnergy *
                   mGaussTurbulentKinematicViscosityDerivative /
                   std::pow(mrData.mTurbulentKinematicViscosity, 2);
    } else {
        return 0.0;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
double KElementDataDerivatives<TDim, TNumNodes>::ShapeDerivative::CalculateSourceTermDerivative() const
{
    return KEpsilonElementData::AdjointUtilities<TDim, TNumNodes>::CalculateProductionShapeDerivative(
        this->mrData.mTurbulentKinematicViscosity, mGaussTurbulentKinematicViscosityDerivative, this->mrData.mSourceTerm,
        this->mrData.mNodalVelocity, this->mrdNdX, this->mrdNdXDerivative);
}

// template instantiations
template class KElementDataDerivatives<2, 3>;
template class KElementDataDerivatives<2, 4>;

template class KElementDataDerivatives<3, 4>;
template class KElementDataDerivatives<3, 8>;

} // namespace KOmegaSSTElementData

} // namespace Kratos
