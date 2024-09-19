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
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/brep_surface.h"
#include "geometries/triangle_2d_3.h"

#include "integration/integration_point_utilities.h"
#include "utilities/quadrature_points_utility.h"

#include <algorithm>

namespace Kratos
{
namespace IgaMappingIntersectionUtilities
{
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    typedef signed long long cInt;

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

    void KRATOS_API(MAPPING_APPLICATION) CreateFEMIgaSurfaceCouplingGeometries(
        ModelPart& rModelPartDomainA,
        ModelPart& rModelPartDomainB,
        ModelPart& rModelPartResult,
        bool origin_is_iga);

    void KRATOS_API(MAPPING_APPLICATION) IgaCreateQuadraturePointsCoupling1DGeometries2D(
        ModelPart& rModelPartCoupling,
        double Tolerance);

    void KRATOS_API(MAPPING_APPLICATION) IgaFEMCreateQuadraturePointsCoupling2DGeometries3D(
        ModelPart& rModelPartCoupling,
        bool origin_is_iga);

    void KRATOS_API(MAPPING_APPLICATION) CreateIgaIgaQuadraturePointsCoupling2DGeometries3D(
        ModelPart& rModelPartCoupling);

    void KRATOS_API(MAPPING_APPLICATION) Triangulation(
        std::vector<CoordinatesArrayType> original_triangle_coordinates,
        GeometryPointerType master_geometry,
        std::vector<std::vector<CoordinatesArrayType>>& new_triangles);

    bool KRATOS_API(MAPPING_APPLICATION) isTriangleIntersectedByKnotLine(
        std::vector<CoordinatesArrayType> triangle_coordinates, 
        double knot_line_position, 
        bool is_vertical);

    void KRATOS_API(MAPPING_APPLICATION) splitTriangle(
        std::vector<std::vector<CoordinatesArrayType>>& new_triangles,
        std::vector<CoordinatesArrayType> triangle_coordinates, 
        double knot_line_position, 
        bool is_vertical);
    
    void KRATOS_API(MAPPING_APPLICATION) intersectionPoint (
        CoordinatesArrayType point_1,
        CoordinatesArrayType point_2, 
        CoordinatesArrayType& intersection_point,
        double knot_line_position, 
        bool is_vertical);

    void KRATOS_API(MAPPING_APPLICATION) categorizePoint(CoordinatesArrayType point,
        std::vector<CoordinatesArrayType>& points_below,
        std::vector<CoordinatesArrayType>& points_above,
        double knot_line_position,
        bool is_vertical
        );
    
    void KRATOS_API(MAPPING_APPLICATION) sortVerticesCounterClockwise (
        std::vector<CoordinatesArrayType>& triangle_vertices,
        std::vector<CoordinatesArrayType>& sorted_triangle_vertices
        );

    void KRATOS_API(MAPPING_APPLICATION) FindIntersectionTriangleSideWithSurfacePatchBisectionMethod (
        GeometryPointerType geom_master,
        CoordinatesArrayType triangle_node_1,
        CoordinatesArrayType triangle_node_2,
        CoordinatesArrayType initial_guess,
        CoordinatesArrayType& intersection_point
        );

    void KRATOS_API(MAPPING_APPLICATION) SubdivideRectangleWithMasterKnotLines(
        std::vector<CoordinatesArrayType> original_rectangle_coordinates,
        GeometryPointerType master_geometry,
        std::vector<std::vector<CoordinatesArrayType>>& new_rectangles);

    void KRATOS_API(MAPPING_APPLICATION) SplitRectangle(
        std::vector<std::vector<CoordinatesArrayType>>& new_rectangles,
        std::vector<CoordinatesArrayType> rectangle_coordinates, 
        double knot_line_position, 
        bool is_vertical);

    bool KRATOS_API(MAPPING_APPLICATION) IsRectangleIntersectedByKnotLine(
        std::vector<CoordinatesArrayType> rectangle_coordinates, 
        double knot_line_position, 
        bool is_vertical);

    void KRATOS_API(MAPPING_APPLICATION) FindInitialGuessNewtonRaphsonProjection(
        CoordinatesArrayType slave_element_node,
        GeometryPointerType master_geometry,
        CoordinatesArrayType& initial_guess);

    bool KRATOS_API(MAPPING_APPLICATION) AreProjectionsOnParameterSpaceBoundary(
        std::vector<CoordinatesArrayType>& points_to_triangulate,
        NurbsSurfaceGeometry<3, PointerVector<Node>>& nurbs_surface);




}  // namespace IgaMappingIntersectionUtilities. 

}  // namespace Kratos.
