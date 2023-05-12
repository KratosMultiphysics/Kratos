//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//                   Peter Wilson
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"

#include "geometries/coupling_geometry.h"

#include "integration/integration_point_utilities.h"
#include "utilities/quadrature_points_utility.h"

namespace Kratos
{
namespace MappingIntersectionUtilities
{
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    typedef Node NodeType;
    typedef typename NodeType::Pointer NodePointerType;
    typedef Geometry<NodeType> GeometryType;
    typedef typename GeometryType::Pointer GeometryPointerType;

    typedef typename GeometryType::GeometriesArrayType GeometriesArrayType;
    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;
    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;

    void KRATOS_API(MAPPING_APPLICATION) FindIntersection1DGeometries2D(
        ModelPart& rModelPartDomainA,
        ModelPart& rModelPartDomainB,
        ModelPart& rModelPartResult,
        double Tolerance = 1e-6);

    void KRATOS_API(MAPPING_APPLICATION) CreateQuadraturePointsCoupling1DGeometries2D(
        ModelPart& rModelPartCoupling,
        double Tolerance);

    bool KRATOS_API(MAPPING_APPLICATION) FindOverlapExtents1DGeometries2D(
        const GeometryType& rMasterLine,
        const GeometryType& rSlaveLine,
        std::vector<array_1d<double, 3 > >& rOverlapExtents,
        const double Tolerance = 1e-6);
}  // namespace MappingIntersectionUtilities.

}  // namespace Kratos.
