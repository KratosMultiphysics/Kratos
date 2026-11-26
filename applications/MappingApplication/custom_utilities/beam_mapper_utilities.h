//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
// Main authors:    E. G. Loera Villeda
// Contributor:     Juan I. Camarotti
//

#pragma once

// System includes

// External includes

// Project includes
#include "projection_utilities.h"

namespace Kratos
{
namespace BeamMapperUtilities
{
    using GeometryType = Geometry<Node>;

    /// Compute Hermitian shape functions and derivatives
    void KRATOS_API(MAPPING_APPLICATION) HermitianShapeFunctionsValues(
        Vector& rHermitianShapeFunctions, 
        Vector& rHermitianShapeFunctionsDer, 
        const array_1d<double, 3>& rCoordinates);

    /// Project a point onto a beam line using Hermitian shape functions
    ProjectionUtilities::PairingIndex KRATOS_API(MAPPING_APPLICATION) ProjectOnLineHermitian(
        const GeometryType& rGeometry,
        const Point& rPointToProject,
        const double LocalCoordTol,
        Vector& rHermitianShapeFunctionValues,
        Vector& rHermitianShapeFunctionValuesDer,
        double& rProjectionDistance,
        Point& rProjectionOfPoint);

}  // namespace BeamMapperUtilities.

}  // namespace Kratos.
