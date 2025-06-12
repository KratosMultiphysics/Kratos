//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#pragma once

// System includes

// External includes

// Project includes
#include "custom_utilities/mapper_utilities.h"
#include "utilities/geometrical_projection_utilities.h"

namespace Kratos
{
namespace ProjectionUtilities
{

enum class PairingIndex
{
    Volume_Inside   = -1,
    Volume_Outside  = -2,
    Surface_Inside  = -3,
    Surface_Outside = -4,
    Line_Inside     = -5,
    Line_Outside    = -6,
    Closest_Point   = -7,
    Unspecified     = -8
};

typedef std::size_t SizeType;
typedef std::size_t IndexType;

typedef Geometry<Node> GeometryType;

PairingIndex KRATOS_API(MAPPING_APPLICATION) ProjectOnLine(const GeometryType& rGeometry,
                           const Point& rPointToProject,
                           const double LocalCoordTol,
                           Vector& rShapeFunctionValues,
                           std::vector<int>& rEquationIds,
                           double& rProjectionDistance,
                           const bool ComputeApproximation=true);

PairingIndex KRATOS_API(MAPPING_APPLICATION) ProjectOnSurface(const GeometryType& rGeometry,
                     const Point& rPointToProject,
                     const double LocalCoordTol,
                     Vector& rShapeFunctionValues,
                     std::vector<int>& rEquationIds,
                     double& rProjectionDistance,
                     const bool ComputeApproximation=true);

PairingIndex KRATOS_API(MAPPING_APPLICATION) ProjectIntoVolume(const GeometryType& rGeometry,
                               const Point& rPointToProject,
                               const double LocalCoordTol,
                               Vector& rShapeFunctionValues,
                               std::vector<int>& rEquationIds,
                               double& rProjectionDistance,
                               const bool ComputeApproximation=true);

bool KRATOS_API(MAPPING_APPLICATION) ComputeProjection(const GeometryType& rGeometry,
                       const Point& rPointToProject,
                       const double LocalCoordTol,
                       Vector& rShapeFunctionValues,
                       std::vector<int>& rEquationIds,
                       double& rProjectionDistance,
                       PairingIndex& rPairingIndex,
                       const bool ComputeApproximation=true);

}  // namespace ProjectionUtilities.

}  // namespace Kratos.
