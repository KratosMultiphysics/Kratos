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

// Application includes
#include "custom_utilities/rans_adjoint_utilities.h"

// Include base h
#include "vms_monolithic_k_based_wall_condition_derivative_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <unsigned int TDim>
array_1d<double, 3> VMSMonolithicKBasedWallConditionDerivativeUtilities<TDim>::VelocityDerivative::CalculateWallVelocityDerivative(
    const IndexType NodeIndex,
    const IndexType DirectionIndex,
    const Vector& rN)
{
    array_1d<double, 3> result = ZeroVector(3);
    result[DirectionIndex] = rN[NodeIndex];
    return result;
}


template <unsigned int TDim>
double VMSMonolithicKBasedWallConditionDerivativeUtilities<TDim>::ShapeDerivative::CalculateWallHeightConditionDerivative(
    const GeometryType& rConditionGeometry,
    const GeometryType& rParentElementGeometry,
    const IndexType NodeIndex,
    const IndexType DirectionIndex,
    const double NormalMagnitude,
    const ArrayD& rUnitNormal)
{
    ArrayD normal_derivative = ZeroVector(3);
    const Vector& temp = row(rConditionGeometry.GetValue(NORMAL_SHAPE_DERIVATIVE),
                             NodeIndex * TDim + DirectionIndex);
    for (IndexType i = 0; i < TDim; ++i) {
        normal_derivative[i] = temp[i];
    }

    const auto& unit_normal_derivative = RansAdjointUtilities::CalculateUnitVectorDerivative(
        NormalMagnitude, rUnitNormal, normal_derivative);

    return RansAdjointUtilities::CalculateWallHeightConditionDerivative(
        rConditionGeometry, rParentElementGeometry, DirectionIndex, rUnitNormal,
        unit_normal_derivative);
}

template <unsigned int TDim>
double VMSMonolithicKBasedWallConditionDerivativeUtilities<TDim>::ShapeDerivative::CalculateWallHeightParentElementDerivative(
    const GeometryType& rConditionGeometry,
    const GeometryType& rParentElementGeometry,
    const IndexType DirectionIndex,
    const ArrayD& rUnitNormal)
{
    return RansAdjointUtilities::CalculateWallHeightParentElementDerivative(
        rConditionGeometry, rParentElementGeometry, DirectionIndex, rUnitNormal,
        ZeroVector(3));
}

// template instantiations
template class VMSMonolithicKBasedWallConditionDerivativeUtilities<2>;
template class VMSMonolithicKBasedWallConditionDerivativeUtilities<3>;

} // namespace Kratos