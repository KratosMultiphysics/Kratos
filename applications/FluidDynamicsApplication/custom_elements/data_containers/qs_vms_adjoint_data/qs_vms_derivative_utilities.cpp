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

template <unsigned int TDim>
void QSVMSDerivativeUtilities<TDim>::CalculateGradient(
    array_1d<double, 3>& rOutput,
    const Variable<double>& rVariable,
    const GeometryType& rGeometry,
    const Matrix& rdNdX)
{
    rOutput.clear();
    for (IndexType a = 0; a < rdNdX.size1(); ++a) {
        const auto& r_node = rGeometry[a];
        for (IndexType i = 0; i < TDim; ++i) {
            rOutput[i] += rdNdX(a, i) * r_node.FastGetSolutionStepValue(rVariable);
        }
    }
}

template <unsigned int TDim>
void QSVMSDerivativeUtilities<TDim>::CalculateGradient(
    BoundedMatrix<double, TDim, TDim>& rOutput,
    const Variable<array_1d<double, 3>>& rVariable,
    const GeometryType& rGeometry,
    const Matrix& rdNdX)
{
    rOutput.clear();
    for (IndexType a = 0; a < rdNdX.size1(); ++a) {
        const auto& r_node = rGeometry[a];
        for (IndexType i = 0; i < TDim; ++i) {
            for (IndexType j = 0; j < TDim; ++j) {
                rOutput(i, j) += rdNdX(a, j) * r_node.FastGetSolutionStepValue(rVariable)[i];
            }
        }
    }
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

template <unsigned int TDim>
QSVMSDerivativeUtilities<TDim>::Derivative::Derivative(
    const IndexType NodeIndex,
    const IndexType DirectionIndex,
    const GeometryType& rGeometry,
    const double W,
    const Vector& rN,
    const Matrix& rdNdX,
    const double WDerivative,
    const double DetJDerivative,
    const Matrix& rdNdXDerivative)
    : mNodeIndex(NodeIndex),
      mDirectionIndex(DirectionIndex),
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

    KRATOS_DEBUG_ERROR_IF(mDirectionIndex >= TDim)
        << "Derivative direction index is not valid. [ mDirectionIndex >= "
           "TDim, mDirectionIndex = "
        << mDirectionIndex << ", TDim = " << TDim << " ].\n";

    KRATOS_CATCH("");
}

template <unsigned int TDim>
template <unsigned int TNumNodes>
array_1d<double, 3> QSVMSDerivativeUtilities<TDim>::VelocityDerivative<TNumNodes>::CalculateEffectiveVelocityDerivative(
    const array_1d<double, 3>& rVelocity) const
{
    array_1d<double, 3> result = ZeroVector(3);
    result[this->mDirectionIndex] = this->mrN[this->mNodeIndex];
    return result;
}

template <unsigned int TDim>
template <unsigned int TNumNodes>
const Variable<double>& QSVMSDerivativeUtilities<TDim>::VelocityDerivative<TNumNodes>::GetDerivativeVariable() const
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

template <unsigned int TDim>
template <unsigned int TNumNodes>
double QSVMSDerivativeUtilities<TDim>::VelocityDerivative<TNumNodes>::CalculateElementLengthDerivative(
    const double ElementLength) const
{
    return 0.0;
}

template <unsigned int TDim>
template <unsigned int TNumNodes>
void QSVMSDerivativeUtilities<TDim>::VelocityDerivative<TNumNodes>::CalculateStrainRateDerivative(
    Vector& rOutput,
    const Matrix& rNodalVelocity) const
{
    QSVMSDerivativeUtilities<TDim>::CalculateStrainRateVelocityDerivative(
        rOutput, this->mNodeIndex, this->mDirectionIndex, this->mrdNdX);
}

template <unsigned int TDim>
template <unsigned int TNumNodes>
array_1d<double, 3> QSVMSDerivativeUtilities<TDim>::PressureDerivative<TNumNodes>::CalculateEffectiveVelocityDerivative(
    const array_1d<double, 3>& rVelocity) const
{
    array_1d<double, 3> result = ZeroVector(3);
    return result;
}

template <unsigned int TDim>
template <unsigned int TNumNodes>
const Variable<double>& QSVMSDerivativeUtilities<TDim>::PressureDerivative<TNumNodes>::GetDerivativeVariable() const
{
    return PRESSURE;
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
template <unsigned int TNumNodes>
array_1d<double, 3> QSVMSDerivativeUtilities<TDim>::ShapeDerivative<TNumNodes>::CalculateEffectiveVelocityDerivative(
    const array_1d<double, 3>& rVelocity) const
{
    array_1d<double, 3> result = ZeroVector(3);
    return result;
}

template <unsigned int TDim>
template <unsigned int TNumNodes>
const Variable<double>& QSVMSDerivativeUtilities<TDim>::ShapeDerivative<TNumNodes>::GetDerivativeVariable() const
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

template <unsigned int TDim>
template <unsigned int TNumNodes>
double QSVMSDerivativeUtilities<TDim>::ShapeDerivative<TNumNodes>::CalculateElementLengthDerivative(
    const double ElementLength) const
{
    return ElementSizeCalculator<TDim, TNumNodes>::MinimumElementSizeDerivative(this->mNodeIndex, this->mDirectionIndex, this->mrGeometry);
}

template <unsigned int TDim>
template <unsigned int TNumNodes>
void QSVMSDerivativeUtilities<TDim>::ShapeDerivative<TNumNodes>::CalculateStrainRateDerivative(
    Vector& rOutput,
    const Matrix& rNodalVelocity) const
{
    CalculateStrainRate(rOutput, rNodalVelocity, this->mrdNdXDerivative);
}

// template instantiations

template class QSVMSDerivativeUtilities<2>::VelocityDerivative<3>;
template class QSVMSDerivativeUtilities<2>::VelocityDerivative<4>;

template class QSVMSDerivativeUtilities<3>::VelocityDerivative<4>;
template class QSVMSDerivativeUtilities<3>::VelocityDerivative<8>;

template class QSVMSDerivativeUtilities<2>::PressureDerivative<3>;
template class QSVMSDerivativeUtilities<2>::PressureDerivative<4>;

template class QSVMSDerivativeUtilities<3>::PressureDerivative<4>;
template class QSVMSDerivativeUtilities<3>::PressureDerivative<8>;

template class QSVMSDerivativeUtilities<2>::ShapeDerivative<3>;
template class QSVMSDerivativeUtilities<2>::ShapeDerivative<4>;

template class QSVMSDerivativeUtilities<3>::ShapeDerivative<4>;
template class QSVMSDerivativeUtilities<3>::ShapeDerivative<8>;

template class QSVMSDerivativeUtilities<2>;
template class QSVMSDerivativeUtilities<3>;



} // namespace Kratos
