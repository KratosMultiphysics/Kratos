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
#include "includes/define.h"
#include "includes/variables.h"
#include "utilities/element_size_calculator.h"

// Application includes

// Include base h
#include "qs_vms_derivative_utilities.h"

namespace Kratos
{
namespace QSVMSDerivativeHelperUtilities
{
template<unsigned int TComponentIndex>
const Variable<double>& GetVelocityVariable();

template<unsigned int TComponentIndex>
const Variable<double>& GetShapeVariable();

template<>
const Variable<double>& GetVelocityVariable<0>()
{
    return VELOCITY_X;
}

template<>
const Variable<double>& GetVelocityVariable<1>()
{
    return VELOCITY_Y;
}

template<>
const Variable<double>& GetVelocityVariable<2>()
{
    return VELOCITY_Z;
}

template<>
const Variable<double>& GetShapeVariable<0>()
{
    return SHAPE_SENSITIVITY_X;
}

template<>
const Variable<double>& GetShapeVariable<1>()
{
    return SHAPE_SENSITIVITY_Y;
}

template<>
const Variable<double>& GetShapeVariable<2>()
{
    return SHAPE_SENSITIVITY_Z;
}
}

template <>
void QSVMSDerivativeUtilities<2>::CalculateStrainRate(
    Vector& rOutput,
    const Matrix& rNodalVelocity,
    const Matrix& rdNdX)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(rOutput.size() != 3)
        << "rOutput is not properly initialized. [ rOutput.size() != 3, "
           "rOutput.size() = "
        << rOutput.size() << " ].\n";

    const IndexType number_of_nodes = rNodalVelocity.size1();

    KRATOS_DEBUG_ERROR_IF(rdNdX.size1() != number_of_nodes)
        << "rdNdX not initialized properly. [ required rdNdX.size1() = " << number_of_nodes
        << ", rdNdX.size1() = " << rdNdX.size1() << " ].\n";

    KRATOS_DEBUG_ERROR_IF(rNodalVelocity.size2() != 2)
        << "rNodalVelocity.size2() != TDim [ TDim = " << 2
        << ", rNodalVelocity.size2() = " << rNodalVelocity.size2() << " ].\n";

    KRATOS_DEBUG_ERROR_IF(rdNdX.size2() != 2)
        << "rdNdX.size2() != TDim [ TDim = " << 2
        << ", rdNdX.size2() = " << rdNdX.size2() << " ].\n";

    noalias(rOutput) = ZeroVector(3);
    for (IndexType i = 0; i < number_of_nodes; ++i) {
        rOutput[0] += rdNdX(i, 0) * rNodalVelocity(i, 0);
        rOutput[1] += rdNdX(i, 1) * rNodalVelocity(i, 1);
        rOutput[2] += rdNdX(i, 0) * rNodalVelocity(i, 1) + rdNdX(i, 1) * rNodalVelocity(i, 0);
    }

    KRATOS_CATCH("");
}

template <>
void QSVMSDerivativeUtilities<3>::CalculateStrainRate(
    Vector& rOutput,
    const Matrix& rNodalVelocity,
    const Matrix& rdNdX)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(rOutput.size() != 6)
        << "rOutput is not properly initialized. [ rOutput.size() != 6, "
           "rOutput.size() = "
        << rOutput.size() << " ].\n";

    const IndexType number_of_nodes = rNodalVelocity.size1();

    KRATOS_DEBUG_ERROR_IF(rdNdX.size1() != number_of_nodes)
        << "rdNdX not initialized properly. [ required rdNdX.size1() = " << number_of_nodes
        << ", rdNdX.size1() = " << rdNdX.size1() << " ].\n";

    KRATOS_DEBUG_ERROR_IF(rNodalVelocity.size2() != 3)
        << "rNodalVelocity.size2() != TDim [ TDim = " << 3
        << ", rNodalVelocity.size2() = " << rNodalVelocity.size2() << " ].\n";

    KRATOS_DEBUG_ERROR_IF(rdNdX.size2() != 3)
        << "rdNdX.size2() != TDim [ TDim = " << 3
        << ", rdNdX.size2() = " << rdNdX.size2() << " ].\n";

    noalias(rOutput) = ZeroVector(6);
    for (IndexType i = 0; i < number_of_nodes; ++i) {
        rOutput[0] += rdNdX(i,0) * rNodalVelocity(i,0);
        rOutput[1] += rdNdX(i,1) * rNodalVelocity(i,1);
        rOutput[2] += rdNdX(i,2) * rNodalVelocity(i,2);
        rOutput[3] += rdNdX(i,0) * rNodalVelocity(i,1) + rdNdX(i,1) * rNodalVelocity(i,0);
        rOutput[4] += rdNdX(i,1) * rNodalVelocity(i,2) + rdNdX(i,2) * rNodalVelocity(i,1);
        rOutput[5] += rdNdX(i,0) * rNodalVelocity(i,2) + rdNdX(i,2) * rNodalVelocity(i,0);
    }

    KRATOS_CATCH("");
}

template <>
void QSVMSDerivativeUtilities<2>::CalculateStrainRateVelocityDerivative(
    Vector& rOutput,
    const IndexType DerivativeNodeIndex,
    const IndexType DerivativeDirectionIndex,
    const Matrix& rdNdX)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(rOutput.size() != 3)
        << "rOutput is not properly initialized. [ rOutput.size() != 3, "
           "rOutput.size() = "
        << rOutput.size() << " ].\n";

    noalias(rOutput) = ZeroVector(3);

    rOutput[DerivativeDirectionIndex] += rdNdX(DerivativeNodeIndex, DerivativeDirectionIndex);
    rOutput[2] += rdNdX(DerivativeNodeIndex, 0) * (DerivativeDirectionIndex == 1);
    rOutput[2] += rdNdX(DerivativeNodeIndex, 1) * (DerivativeDirectionIndex == 0);

    KRATOS_CATCH("");
}

template <>
void QSVMSDerivativeUtilities<3>::CalculateStrainRateVelocityDerivative(
    Vector& rOutput,
    const IndexType DerivativeNodeIndex,
    const IndexType DerivativeDirectionIndex,
    const Matrix& rdNdX)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(rOutput.size() != 6)
        << "rOutput is not properly initialized. [ rOutput.size() != 6, "
           "rOutput.size() = "
        << rOutput.size() << " ].\n";

    noalias(rOutput) = ZeroVector(6);
    rOutput[DerivativeDirectionIndex] += rdNdX(DerivativeNodeIndex, DerivativeDirectionIndex);

    rOutput[3] += rdNdX(DerivativeNodeIndex, 0) * (DerivativeDirectionIndex == 1);
    rOutput[3] += rdNdX(DerivativeNodeIndex, 1) * (DerivativeDirectionIndex == 0);

    rOutput[4] += rdNdX(DerivativeNodeIndex, 1) * (DerivativeDirectionIndex == 2);
    rOutput[4] += rdNdX(DerivativeNodeIndex, 2) * (DerivativeDirectionIndex == 1);

    rOutput[5] += rdNdX(DerivativeNodeIndex, 0) * (DerivativeDirectionIndex == 2);
    rOutput[5] += rdNdX(DerivativeNodeIndex, 2) * (DerivativeDirectionIndex == 0);

    KRATOS_CATCH("");
}

template<>
const std::array<const Variable<double>*, 3> QSVMSDerivativeUtilities<2>::GetStrainRateVariables()
{
    return {&STRAIN_RATE_2D_XX, &STRAIN_RATE_2D_YY, &STRAIN_RATE_2D_XY};
}

template<>
const std::array<const Variable<double>*, 6> QSVMSDerivativeUtilities<3>::GetStrainRateVariables()
{
    return {&STRAIN_RATE_3D_XX, &STRAIN_RATE_3D_YY, &STRAIN_RATE_3D_ZZ, &STRAIN_RATE_3D_XY, &STRAIN_RATE_3D_YZ, &STRAIN_RATE_3D_XZ};
}

template <unsigned int TDim>
template <unsigned int TComponentIndex>
QSVMSDerivativeUtilities<TDim>::Derivative<TComponentIndex>::Derivative(
    const IndexType NodeIndex,
    const GeometryType& rGeometry,
    const double W,
    const Vector& rN,
    const Matrix& rdNdX,
    const double WDerivative,
    const double DetJDerivative,
    const Matrix& rdNdXDerivative)
    : mNodeIndex(NodeIndex),
      mrGeometry(rGeometry),
      mW(W),
      mrN(rN),
      mrdNdX(rdNdX),
      mWDerivative(WDerivative),
      mDetJDerivative(DetJDerivative),
      mrdNdXDerivative(rdNdXDerivative)
{
    KRATOS_TRY

    const IndexType number_of_nodes = rGeometry.PointsNumber();

    KRATOS_DEBUG_ERROR_IF(mrN.size() != number_of_nodes)
        << "mrN vector is not initialized properly. [ mrN.size() != "
           "number_of_nodes, mrN.size() = "
        << mrN.size() << ", number_of_nodes = " << number_of_nodes << " ].\n";

    KRATOS_DEBUG_ERROR_IF(mrdNdX.size1() != number_of_nodes)
        << "mrdNdX matrix is not initialized properly. [ mrdNdX.size1() != "
           "number_of_nodes, mrdNdX.size1() = "
        << mrdNdX.size1() << ", number_of_nodes = " << number_of_nodes << " ].\n";

    KRATOS_DEBUG_ERROR_IF(mrdNdX.size2() != TDim)
        << "mrdNdX matrix is not initialized properly. [ mrdNdX.size2() != "
           "TDim, mrdNdX.size2() = "
        << mrdNdX.size2() << ", TDim = " << TDim << " ].\n";

    KRATOS_DEBUG_ERROR_IF(mrdNdXDerivative.size1() != number_of_nodes)
        << "mrdNdXDerivative matrix is not initialized properly. [ mrdNdXDerivative.size1() != "
           "number_of_nodes, mrdNdXDerivative.size1() = "
        << mrdNdXDerivative.size1() << ", number_of_nodes = " << number_of_nodes << " ].\n";

    KRATOS_DEBUG_ERROR_IF(mrdNdXDerivative.size2() != TDim)
        << "mrdNdXDerivative matrix is not initialized properly. [ mrdNdXDerivative.size2() != "
           "TDim, mrdNdXDerivative.size2() = "
        << mrdNdXDerivative.size2() << ", TDim = " << TDim << " ].\n";

    KRATOS_DEBUG_ERROR_IF(mNodeIndex >= number_of_nodes)
        << "Derivative node index is not valid. [ mNodeIndex >= "
           "number_of_nodes, mNodeIndex = "
        << mNodeIndex << ", number_of_nodes = " << number_of_nodes << " ].\n";

    KRATOS_CATCH("");
}

template <unsigned int TDim>
template <unsigned int TNumNodes, unsigned int TComponentIndex>
const Variable<double>& QSVMSDerivativeUtilities<TDim>::VelocityDerivative<TNumNodes, TComponentIndex>::GetDerivativeVariable() const
{
    static_assert(TDim > TComponentIndex);
    return QSVMSDerivativeHelperUtilities::GetVelocityVariable<TComponentIndex>();
}

template <unsigned int TDim>
template <unsigned int TNumNodes, unsigned int TComponentIndex>
array_1d<double, TDim> QSVMSDerivativeUtilities<TDim>::VelocityDerivative<TNumNodes, TComponentIndex>::CalculateEffectiveVelocityDerivative(
    const array_1d<double, TDim>& rVelocity) const
{
    array_1d<double, TDim> result = ZeroVector(TDim);
    result[TComponentIndex] = this->mrN[this->mNodeIndex];
    return result;
}

template <unsigned int TDim>
template <unsigned int TNumNodes, unsigned int TComponentIndex>
double QSVMSDerivativeUtilities<TDim>::VelocityDerivative<TNumNodes, TComponentIndex>::CalculateElementLengthDerivative(
    const double ElementLength) const
{
    return 0.0;
}

template <unsigned int TDim>
template <unsigned int TNumNodes, unsigned int TComponentIndex>
void QSVMSDerivativeUtilities<TDim>::VelocityDerivative<TNumNodes, TComponentIndex>::CalculateStrainRateDerivative(
    Vector& rOutput,
    const Matrix& rNodalVelocity) const
{
    QSVMSDerivativeUtilities<TDim>::CalculateStrainRateVelocityDerivative(
        rOutput, this->mNodeIndex, TComponentIndex, this->mrdNdX);
}

template <unsigned int TDim>
template <unsigned int TNumNodes>
array_1d<double, TDim> QSVMSDerivativeUtilities<TDim>::PressureDerivative<TNumNodes>::CalculateEffectiveVelocityDerivative(
    const array_1d<double, TDim>& rVelocity) const
{
    array_1d<double, TDim> result = ZeroVector(TDim);
    return result;
}

template <unsigned int TDim>
template <unsigned int TNumNodes>
double QSVMSDerivativeUtilities<TDim>::PressureDerivative<TNumNodes>::CalculateElementLengthDerivative(
    const double ElementLength) const
{
    return 0.0;
}

template <unsigned int TDim>
template <unsigned int TNumNodes>
void QSVMSDerivativeUtilities<TDim>::PressureDerivative<TNumNodes>::CalculateStrainRateDerivative(
    Vector& rOutput,
    const Matrix& rNodalVelocity) const
{
    rOutput.clear();
}

template <unsigned int TDim>
template <unsigned int TNumNodes, unsigned int TComponentIndex>
const Variable<double>& QSVMSDerivativeUtilities<TDim>::ShapeDerivative<TNumNodes, TComponentIndex>::GetDerivativeVariable() const
{
    static_assert(TDim > TComponentIndex);
    return QSVMSDerivativeHelperUtilities::GetShapeVariable<TComponentIndex>();
}

template <unsigned int TDim>
template <unsigned int TNumNodes, unsigned int TComponentIndex>
array_1d<double, TDim> QSVMSDerivativeUtilities<TDim>::ShapeDerivative<TNumNodes, TComponentIndex>::CalculateEffectiveVelocityDerivative(
    const array_1d<double, TDim>& rVelocity) const
{
    array_1d<double, TDim> result = ZeroVector(TDim);
    return result;
}

template <unsigned int TDim>
template <unsigned int TNumNodes, unsigned int TComponentIndex>
double QSVMSDerivativeUtilities<TDim>::ShapeDerivative<TNumNodes, TComponentIndex>::CalculateElementLengthDerivative(
    const double ElementLength) const
{
    return ElementSizeCalculator<TDim, TNumNodes>::MinimumElementSizeDerivative(this->mNodeIndex, TComponentIndex, this->mrGeometry);
}

template <unsigned int TDim>
template <unsigned int TNumNodes, unsigned int TComponentIndex>
void QSVMSDerivativeUtilities<TDim>::ShapeDerivative<TNumNodes, TComponentIndex>::CalculateStrainRateDerivative(
    Vector& rOutput,
    const Matrix& rNodalVelocity) const
{
    CalculateStrainRate(rOutput, rNodalVelocity, this->mrdNdXDerivative);
}

// template instantiations
template class QSVMSDerivativeUtilities<2>::Derivative<0>;
template class QSVMSDerivativeUtilities<2>::Derivative<1>;

template class QSVMSDerivativeUtilities<3>::Derivative<0>;
template class QSVMSDerivativeUtilities<3>::Derivative<1>;
template class QSVMSDerivativeUtilities<3>::Derivative<2>;

template class QSVMSDerivativeUtilities<2>::VelocityDerivative<3, 0>;
template class QSVMSDerivativeUtilities<2>::VelocityDerivative<3, 1>;

template class QSVMSDerivativeUtilities<2>::VelocityDerivative<4, 0>;
template class QSVMSDerivativeUtilities<2>::VelocityDerivative<4, 1>;

template class QSVMSDerivativeUtilities<3>::VelocityDerivative<4, 0>;
template class QSVMSDerivativeUtilities<3>::VelocityDerivative<4, 1>;
template class QSVMSDerivativeUtilities<3>::VelocityDerivative<4, 2>;

template class QSVMSDerivativeUtilities<3>::VelocityDerivative<8, 0>;
template class QSVMSDerivativeUtilities<3>::VelocityDerivative<8, 1>;
template class QSVMSDerivativeUtilities<3>::VelocityDerivative<8, 2>;

template class QSVMSDerivativeUtilities<2>::PressureDerivative<3>;
template class QSVMSDerivativeUtilities<2>::PressureDerivative<4>;

template class QSVMSDerivativeUtilities<3>::PressureDerivative<4>;
template class QSVMSDerivativeUtilities<3>::PressureDerivative<8>;

template class QSVMSDerivativeUtilities<2>::ShapeDerivative<3, 0>;
template class QSVMSDerivativeUtilities<2>::ShapeDerivative<3, 1>;

template class QSVMSDerivativeUtilities<2>::ShapeDerivative<4, 0>;
template class QSVMSDerivativeUtilities<2>::ShapeDerivative<4, 1>;

template class QSVMSDerivativeUtilities<3>::ShapeDerivative<4, 0>;
template class QSVMSDerivativeUtilities<3>::ShapeDerivative<4, 1>;
template class QSVMSDerivativeUtilities<3>::ShapeDerivative<4, 2>;

template class QSVMSDerivativeUtilities<3>::ShapeDerivative<8, 0>;
template class QSVMSDerivativeUtilities<3>::ShapeDerivative<8, 1>;
template class QSVMSDerivativeUtilities<3>::ShapeDerivative<8, 2>;

template class QSVMSDerivativeUtilities<2>;
template class QSVMSDerivativeUtilities<3>;

} // namespace Kratos
