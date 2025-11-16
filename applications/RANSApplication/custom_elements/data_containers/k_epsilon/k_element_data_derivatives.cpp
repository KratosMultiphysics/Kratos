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
#include "custom_utilities/fluid_calculation_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "k_element_data_derivatives.h"

namespace Kratos
{

namespace KEpsilonElementData
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

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENCE_RANS_C_MU))
        << "TURBULENCE_RANS_C_MU is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENT_KINETIC_ENERGY_SIGMA))
        << "TURBULENT_KINETIC_ENERGY_SIGMA is not found in process info.\n";

    KRATOS_ERROR_IF_NOT(r_properties.Has(DYNAMIC_VISCOSITY))
        << "DYNAMIC_VISCOSITY is not found in element properties [ Element.Id() = "
        << rElement.Id() << ", Properties.Id() = " << r_properties.Id() << " ].\n";
    KRATOS_ERROR_IF_NOT(r_properties.Has(DENSITY))
        << "DENSITY is not found in element properties [ Element.Id() = "
        << rElement.Id() << ", Properties.Id() = " << r_properties.Id() << " ].\n";

    for (int i_node = 0; i_node < number_of_nodes; ++i_node) {
        const auto& r_node = r_geometry[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_AUXILIARY_VARIABLE_1, r_node);

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
    using namespace KEpsilonElementData;

    const auto& r_geometry = this->GetGeometry();

    FluidCalculationUtilities::EvaluateInPoint(
        r_geometry, rN, Step,
        std::tie(mTurbulentKineticEnergy, TURBULENT_KINETIC_ENERGY),
        std::tie(mTurbulentEnergyDissipationRate, TURBULENT_ENERGY_DISSIPATION_RATE),
        std::tie(mEffectiveVelocity, VELOCITY));

    FluidCalculationUtilities::EvaluateGradientInPoint(
        r_geometry, rdNdX, Step,
        std::tie(mVelocityGradient, VELOCITY));

    mVelocityDivergence = RansCalculationUtilities::CalculateMatrixTrace<TDim>(mVelocityGradient);

    mTurbulentKinematicViscosity = mCmu * std::pow(mTurbulentKineticEnergy, 2) / mTurbulentEnergyDissipationRate;

    mGamma = CalculateGamma(mCmu, mTurbulentKineticEnergy, mTurbulentKinematicViscosity);

    for (IndexType i = 0; i < TNumNodes; ++i) {
        const auto& velocity = r_geometry[i].FastGetSolutionStepValue(VELOCITY, Step);
        for (IndexType j = 0; j < TDim; ++j) {
            mNodalVelocity(i, j) = velocity[j];
        }
    }

    mEffectiveKinematicViscosity = mKinematicViscosity + mTurbulentKinematicViscosity * mInvTkeSigma;
    mReactionTerm = std::max(mGamma + (2.0 / 3.0) * mVelocityDivergence, 0.0);
    mSourceTerm = KEpsilonElementData::CalculateProductionTerm<TDim>(mVelocityGradient, mTurbulentKinematicViscosity);
}

/***************************************************************/
/********************* Velocity Derivative *********************/
/***************************************************************/

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
    return 0.0;
}

template <unsigned int TDim, unsigned int TNumNodes>
double KElementDataDerivatives<TDim, TNumNodes>::UDerivative::CalculateReactionTermDerivative() const
{
    if (mrData.mReactionTerm > 0.0) {
        return this->mrdNdX(this->mNodeIndex, this->mDirectionIndex) * (2.0 / 3.0);
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
        0.0,
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
    mGaussTurbulentKinematicViscosityDerivative =
        mrData.mCmu * 2.0 * mrData.mTurbulentKineticEnergy /
        mrData.mTurbulentEnergyDissipationRate * this->mrN[this->mNodeIndex];
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
    return mGaussTurbulentKinematicViscosityDerivative * mrData.mInvTkeSigma;
}

template <unsigned int TDim, unsigned int TNumNodes>
double KElementDataDerivatives<TDim, TNumNodes>::KDerivative::CalculateReactionTermDerivative() const
{
    if (mrData.mReactionTerm > 0.0) {
        return KEpsilonElementData::AdjointUtilities<TDim, TNumNodes>::CalculateGammaKDerivative(
            this->mNodeIndex, mrData.mCmu, mrData.mGamma, mrData.mTurbulentKinematicViscosity,
            mGaussTurbulentKinematicViscosityDerivative, this->mrN);
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
/********************** Epsilon Derivative *********************/
/***************************************************************/

template <unsigned int TDim, unsigned int TNumNodes>
KElementDataDerivatives<TDim, TNumNodes>::EpsilonDerivative::EpsilonDerivative(
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
    mGaussTurbulentKinematicViscosityDerivative =
        -1.0 * mrData.mCmu *
        std::pow(mrData.mTurbulentKineticEnergy / mrData.mTurbulentEnergyDissipationRate, 2) *
        this->mrN[this->mNodeIndex];
}

template <unsigned int TDim, unsigned int TNumNodes>
const Variable<double>& KElementDataDerivatives<TDim, TNumNodes>::EpsilonDerivative::GetDerivativeVariable() const
{
    return TURBULENT_ENERGY_DISSIPATION_RATE;
}

template <unsigned int TDim, unsigned int TNumNodes>
array_1d<double, TDim> KElementDataDerivatives<TDim, TNumNodes>::EpsilonDerivative::CalculateEffectiveVelocityDerivative() const
{
    return ZeroVector(TDim);
}

template <unsigned int TDim, unsigned int TNumNodes>
double KElementDataDerivatives<TDim, TNumNodes>::EpsilonDerivative::CalculateEffectiveKinematicViscosityDerivative() const
{
    return mGaussTurbulentKinematicViscosityDerivative * mrData.mInvTkeSigma;
}

template <unsigned int TDim, unsigned int TNumNodes>
double KElementDataDerivatives<TDim, TNumNodes>::EpsilonDerivative::CalculateReactionTermDerivative() const
{
    if (mrData.mReactionTerm > 0.0) {
        return KEpsilonElementData::AdjointUtilities<TDim, TNumNodes>::CalculateGammaEpsilonDerivative(
            mrData.mGamma, mrData.mTurbulentKinematicViscosity,
            mGaussTurbulentKinematicViscosityDerivative);
    } else {
        return 0.0;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
double KElementDataDerivatives<TDim, TNumNodes>::EpsilonDerivative::CalculateSourceTermDerivative() const
{
    return KEpsilonElementData::AdjointUtilities<TDim, TNumNodes>::CalculateProductionScalarDerivative(
        mrData.mTurbulentKinematicViscosity, mrData.mSourceTerm,
        mGaussTurbulentKinematicViscosityDerivative);
}

/***************************************************************/
/*********************** Shape Derivative **********************/
/***************************************************************/

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
    return 0.0;
}

template <unsigned int TDim, unsigned int TNumNodes>
double KElementDataDerivatives<TDim, TNumNodes>::ShapeDerivative::CalculateReactionTermDerivative() const
{
    if (mrData.mReactionTerm > 0.0) {
        MatrixDD gradient_derivative;
        FluidCalculationUtilities::EvaluateGradientInPoint(
            this->mrData.GetGeometry(), this->mrdNdXDerivative,
            std::tie(gradient_derivative, VELOCITY));
        return (2.0 / 3.0) * RansCalculationUtilities::CalculateMatrixTrace<TDim>(gradient_derivative);
    } else {
        return 0.0;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
double KElementDataDerivatives<TDim, TNumNodes>::ShapeDerivative::CalculateSourceTermDerivative() const
{
    return KEpsilonElementData::AdjointUtilities<TDim, TNumNodes>::CalculateProductionShapeDerivative(
        this->mrData.mTurbulentKinematicViscosity, 0.0, this->mrData.mSourceTerm,
        this->mrData.mNodalVelocity, this->mrdNdX, this->mrdNdXDerivative);
}

// template instantiations
template class KElementDataDerivatives<2, 3>;
template class KElementDataDerivatives<2, 4>;

template class KElementDataDerivatives<3, 4>;
template class KElementDataDerivatives<3, 8>;

} // namespace KEpsilonElementData

} // namespace Kratos
