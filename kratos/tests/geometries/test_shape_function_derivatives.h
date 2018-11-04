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
void TestShapeFunctionsLocalGradients(Geometry<Node<3>> const& rGeom);

void TestShapeFunctionsLocalGradients(Geometry<Node<3>> const& rGeom,
                                      GeometryData::IntegrationMethod ThisMethod);

void TestShapeFunctionsLocalGradients_IntegrationPointIndex(Geometry<Node<3>> const& rGeom);

void TestShapeFunctionsLocalGradients_IntegrationPointIndex(
    Geometry<Node<3>> const& rGeom, GeometryData::IntegrationMethod ThisMethod);

void TestShapeFunctionsLocalGradients_IntegrationPointIndex_ShapeFunctionIndex(
    Geometry<Node<3>> const& rGeom, GeometryData::IntegrationMethod ThisMethod);

void TestShapeFunctionsLocalGradients_Point(Geometry<Node<3>> const& rGeom);

void TestAllShapeFunctionsLocalGradients(Geometry<Node<3>> const& rGeom);

void TestShapeFunctionsLocalGradient(Geometry<Node<3>> const& rGeom,
                                     Geometry<Node<3>>::IntegrationPointType Point,
                                     Matrix const& rLocalGradient);
}
}

#endif // KRATOS_TEST_SHAPE_FUNCTION_DERIVATIVES_H_INCLUDED defined
