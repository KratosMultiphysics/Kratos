//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Juan Ignacio Camarotti
//                   Andrea Gorgi
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"

#include "geometries/coupling_geometry.h"
#include "geometries/nurbs_curve_on_surface_geometry.h"

#include "integration/integration_point_utilities.h"
#include "utilities/quadrature_points_utility.h"

namespace Kratos
{
namespace IgaMappingIntersectionUtilities
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

    void KRATOS_API(MAPPING_APPLICATION) IgaCreateBrepCurveOnSurfaceCouplingGeometries(
        ModelPart& rModelPartDomainA,
        ModelPart& rModelPartDomainB,
        ModelPart& rModelPartResult,
        double Tolerance = 1e-6);

    void KRATOS_API(MAPPING_APPLICATION) IgaCreateQuadraturePointsCoupling1DGeometries2D(
        ModelPart& rModelPartCoupling,
        double Tolerance);

}  // namespace IgaMappingIntersectionUtilities. 

}  // namespace Kratos.
