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

DivideTriangle3D3::DivideTriangle3D3(const GeometryType& rInputGeometry, const Vector& rNodalDistances) :
    DivideTriangle2D3(rInputGeometry, rNodalDistances) {};

DivideTriangle3D3::~DivideTriangle3D3() {};

DivideTriangle3D3::IndexedPointGeometryPointerType DivideTriangle3D3::GenerateAuxiliaryPartitionTriangle(
    const int I0,
    const int I1,
    const int I2)
{
    return Kratos::make_shared<IndexedPointTriangleType>(
        mAuxPointsContainer(I0),
        mAuxPointsContainer(I1),
        mAuxPointsContainer(I2));
}

DivideTriangle3D3::IndexedPointGeometryPointerType DivideTriangle3D3::GenerateIntersectionLine(
    const int I0,
    const int I1)
{
    return Kratos::make_shared<IndexedPointLineType>(
        mAuxPointsContainer(I0),
        mAuxPointsContainer(I1));
}

};
