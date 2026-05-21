//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//   License:        BSD License
//                   Kratos default license: kratos/license.txt
//
//   Main authors:   Pablo Becker
//

// System includes

// External includes

// Project includes

#include "utilities/divide_triangle_3d_3.h"

namespace Kratos
{

template<class TPointType>
DivideTriangle3D3<TPointType>::DivideTriangle3D3(const GeometryType& rInputGeometry, const Vector& rNodalDistances) :
    DivideTriangle2D3<TPointType>(rInputGeometry, rNodalDistances) {};

template<class TPointType>
DivideTriangle3D3<TPointType>::~DivideTriangle3D3() {};

template<class TPointType>
typename DivideTriangle3D3<TPointType>::IndexedPointGeometryPointerType DivideTriangle3D3<TPointType>::GenerateAuxiliaryPartitionTriangle(
    const int I0,
    const int I1,
    const int I2)
{
    return Kratos::make_shared<IndexedPointTriangleType>(
        this->mAuxPointsContainer(I0),
        this->mAuxPointsContainer(I1),
        this->mAuxPointsContainer(I2));
}

template<class TPointType>
typename DivideTriangle3D3<TPointType>::IndexedPointGeometryPointerType DivideTriangle3D3<TPointType>::GenerateIntersectionLine(
    const int I0,
    const int I1)
{
    return Kratos::make_shared<IndexedPointLineType>(
        this->mAuxPointsContainer(I0),
        this->mAuxPointsContainer(I1));
}

template class DivideTriangle3D3<Node>;
template class DivideTriangle3D3<IndexedPoint>;

};
