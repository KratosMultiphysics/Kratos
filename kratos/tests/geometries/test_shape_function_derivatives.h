//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    
//
//

#if !defined(KRATOS_TEST_SHAPE_FUNCTION_DERIVATIVES_H_INCLUDED)
#define KRATOS_TEST_SHAPE_FUNCTION_DERIVATIVES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/ublas_interface.h"

namespace Kratos {
namespace Testing {
///@name Type Definitions
///@{

typedef Node<3> NodeType;
typedef Geometry<NodeType> GeometryType;

///@}

void TestShapeFunctionsLocalGradients(GeometryType const& rGeom);

void TestShapeFunctionsLocalGradients(GeometryType const& rGeom,
                                      GeometryData::IntegrationMethod ThisMethod);

void TestShapeFunctionsLocalGradients_IntegrationPointIndex(GeometryType const& rGeom);

void TestShapeFunctionsLocalGradients_IntegrationPointIndex(
    GeometryType const& rGeom, GeometryData::IntegrationMethod ThisMethod);

void TestShapeFunctionsLocalGradients_IntegrationPointIndex_ShapeFunctionIndex(
    GeometryType const& rGeom, GeometryData::IntegrationMethod ThisMethod);

void TestShapeFunctionsLocalGradients_Point(GeometryType const& rGeom);

void TestAllShapeFunctionsLocalGradients(GeometryType const& rGeom);

void TestShapeFunctionsLocalGradient(GeometryType const& rGeom,
                                     GeometryType::IntegrationPointType Point,
                                     Matrix const& rLocalGradient);
}
}

#endif // KRATOS_TEST_SHAPE_FUNCTION_DERIVATIVES_H_INCLUDED defined
