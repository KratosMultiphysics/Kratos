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
namespace IgaMappingIntersectionUtilities
{
    using SizeType = std::size_t;
    using IndexType = std::size_t;

    using NodeType = Node;
    using GeometryType = Geometry<NodeType>;

    using GeometriesArrayType = typename GeometryType::GeometriesArrayType;
    using CoordinatesArrayType = typename GeometryType::CoordinatesArrayType;
    using IntegrationPointsArrayType = typename GeometryType::IntegrationPointsArrayType;

    // This function creates coupling geometries between the IGA and FEM interface. This coupling geometry is composed of a brep_curve_on_surface (IGA side) and 
    // a nurbs curve (FEM side)
    void KRATOS_API(MAPPING_APPLICATION) CreateIgaFEMCouplingGeometries(
        ModelPart& rModelPartDomainA,
        ModelPart& rModelPartDomainB,
        const bool& rIsOriginIga,
        ModelPart& rModelPartResult, 
        double Tolerance = 1.0e-6);
    
    // This function creates quadrature points along the coupling interface (intersection of both domains)
    void KRATOS_API(MAPPING_APPLICATION) CreateIgaFEMQuadraturePointsCouplingInterface(
        ModelPart& rModelPartCoupling,
        double Tolerance);

}  // namespace IgaMappingIntersectionUtilities.

}  // namespace Kratos.
