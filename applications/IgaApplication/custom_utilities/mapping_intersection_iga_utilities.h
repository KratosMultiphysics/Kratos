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
#include "includes/properties.h"

#include "geometries/coupling_geometry.h"
#include "geometries/brep_surface.h"


#include "integration/integration_point_utilities.h"
#include "utilities/quadrature_points_utility.h"

#include "integration/integration_info.h"
#include "integration/triangle_gauss_legendre_integration_points.h"

namespace Kratos
{
namespace MappingIntersectionIgaUtilities
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

    typedef Point EmbeddedNodeType;
    typedef PointerVector<NodeType> ContainerNodeType;
    typedef PointerVector<EmbeddedNodeType> ContainerEmbeddedNodeType;
    typedef BrepSurface<ContainerNodeType, false, ContainerEmbeddedNodeType> BrepSurfaceType;


    void KRATOS_API(MAPPING_APPLICATION) FindIntersection1DGeometries2D(
        ModelPart& rModelPartDomainA,
        ModelPart& rModelPartDomainB,
        ModelPart& rModelPartResult,
        double Tolerance = 1e-6);

    void KRATOS_API(MAPPING_APPLICATION) CreateQuadraturePointsCoupling1DGeometries2D(
        ModelPart& rModelPartCoupling,
        double Tolerance);

    void KRATOS_API(MAPPING_APPLICATION) FindIntersection2DGeometries3D(
        ModelPart& rModelPartDomainA,
        ModelPart& rModelPartDomainB,
        ModelPart& rModelPartResult,
        double Tolerance = 1e-6);

    void KRATOS_API(MAPPING_APPLICATION) CreateQuadraturePointsCoupling2DGeometries3D(
        ModelPart& rModelPartCoupling,
        double Tolerance);
}  // namespace MappingIntersectionUtilities.

}  // namespace Kratos.
