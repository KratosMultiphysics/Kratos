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

#if !defined(KRATOS_CROSS_CHECK_SHAPE_FUNCTIONS_VALUES_H_INCLUDED)
#define KRATOS_CROSS_CHECK_SHAPE_FUNCTIONS_VALUES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/ublas_interface.h"

namespace Kratos {
namespace Testing {
void CrossCheckShapeFunctionsValues(Geometry<Node<3>> const& rGeom);
}
}

#endif // KRATOS_CROSS_CHECK_SHAPE_FUNCTIONS_VALUES_H_INCLUDED defined
