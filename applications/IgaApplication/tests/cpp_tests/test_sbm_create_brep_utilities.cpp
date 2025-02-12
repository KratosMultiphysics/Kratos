//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolo' Antonelli
//  Main authors:    Andrea Gorgi

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "custom_utilities/create_breps_sbm_utilities.h"
#include "includes/kratos_parameters.h"
#include "custom_modelers/nurbs_geometry_modeler.h"

namespace Kratos::Testing
{
    using NurbsSurfaceGeometryType = NurbsSurfaceGeometry<3, PointerVector<Node>>;
    using NurbsSurfaceGeometryPointerType = typename NurbsSurfaceGeometryType::Pointer;

// Tests the SnakeSbmUtilities with a square outer geometry
KRATOS_TEST_CASE_IN_SUITE(TestCreateBrepsSbmUtilitiesOuter, KratosIgaFastSuite)
{
    using IndexType = std::size_t;
    using GeometryType = Geometry<Node>;
    using ContainerNodeType = PointerVector<Node>;
    using ContainerEmbeddedNodeType = PointerVector<Point>;
    using BrepCurveOnSurfaceType = BrepCurveOnSurface<ContainerNodeType, true, ContainerEmbeddedNodeType>;
    using CoordinatesArrayType = GeometryType::CoordinatesArrayType;
    
    // Initialization
    std::size_t EchoLevel = 0;
    Model model;
    ModelPart& iga_model_part = model.CreateModelPart("IgaModelPart");
    ModelPart& surrogate_sub_model_part_outer = iga_model_part.CreateSubModelPart("surrogate_outer");
    ModelPart& surrogate_sub_model_part_inner = iga_model_part.CreateSubModelPart("surrogate_inner");
    
    surrogate_sub_model_part_outer.CreateNewProperties(0);
    
    surrogate_sub_model_part_outer.CreateNewNode(1, -1.0, -1.0, 0.0);
    surrogate_sub_model_part_outer.CreateNewNode(2, -1.0, -0.9, 0.0);
    surrogate_sub_model_part_outer.CreateNewNode(3, -1.0, -0.8, 0.0);

    surrogate_sub_model_part_outer.CreateNewNode(4, -0.9, -0.8, 0.0);
    surrogate_sub_model_part_outer.CreateNewNode(5, -0.8, -0.8, 0.0);
    
    surrogate_sub_model_part_outer.CreateNewNode(6, -0.8, -0.9, 0.0);
    surrogate_sub_model_part_outer.CreateNewNode(7, -0.8, -1.0, 0.0);

    surrogate_sub_model_part_outer.CreateNewNode(8, -0.9, -1.0, 0.0);
    
    
    Properties::Pointer p_prop = surrogate_sub_model_part_outer.pGetProperties(0);
    surrogate_sub_model_part_outer.CreateNewCondition("LineCondition2D2N", 1, {{1, 2}}, p_prop);
    surrogate_sub_model_part_outer.CreateNewCondition("LineCondition2D2N", 2, {{2, 3}}, p_prop);
    surrogate_sub_model_part_outer.CreateNewCondition("LineCondition2D2N", 3, {{3, 4}}, p_prop);
    surrogate_sub_model_part_outer.CreateNewCondition("LineCondition2D2N", 4, {{4, 5}}, p_prop);
    surrogate_sub_model_part_outer.CreateNewCondition("LineCondition2D2N", 5, {{5, 6}}, p_prop);
    surrogate_sub_model_part_outer.CreateNewCondition("LineCondition2D2N", 6, {{6, 7}}, p_prop);
    surrogate_sub_model_part_outer.CreateNewCondition("LineCondition2D2N", 7, {{7, 8}}, p_prop);

    Kratos::Parameters parameters(R"(
        {
            "echo_level": 0,
            "model_part_name": "IgaModelPart",
            "geometry_name"  : "nurbs_surface",
            "lower_point_xyz": [-1, -1, 0],
            "upper_point_xyz": [1, 1, 0],
            "lower_point_uvw": [-1, -1, 0],
            "upper_point_uvw": [1, 1, 0],
            "polynomial_order":     [2, 2],
            "number_of_knot_spans": [20, 20]
        }
    )");

    const Point& A_uvw = Point(-1, -1, 0);
    const Point& B_uvw = Point(1, 1, 0);

    NurbsGeometryModeler nurbs_modeler(model, parameters);
    nurbs_modeler.SetupGeometryModel();
    
    NurbsSurfaceGeometryPointerType surface = std::dynamic_pointer_cast<NurbsSurfaceGeometryType>(iga_model_part.pGetGeometry("nurbs_surface"));
    
    // Create the breps for the outer sbm boundary
    CreateBrepsSbmUtilities<Node, Point> CreateBrepsSbmUtilities(EchoLevel);
    CreateBrepsSbmUtilities.CreateSurrogateBoundary(surface, iga_model_part, surrogate_sub_model_part_inner, surrogate_sub_model_part_outer, A_uvw, B_uvw);
    
    const double tolerance = 1e-12;

    for (IndexType current_id = 2; current_id < iga_model_part.NumberOfGeometries(); current_id++) {

        BrepCurveOnSurfaceType::Pointer brep_curve = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(iga_model_part.pGetGeometry(current_id));
        auto ppp = brep_curve->DomainInterval();
        CoordinatesArrayType vertex_1_on_physical;
        CoordinatesArrayType vertex_2_on_physical;
        CoordinatesArrayType local_coord_vertex_1 = ZeroVector(3); 
        CoordinatesArrayType local_coord_vertex_2 = ZeroVector(3);

        local_coord_vertex_1[0] = ppp.GetT0(); local_coord_vertex_2[0] = ppp.GetT1();

        brep_curve->GlobalCoordinates(vertex_1_on_physical, local_coord_vertex_1);
        brep_curve->GlobalCoordinates(vertex_2_on_physical, local_coord_vertex_2);

        KRATOS_EXPECT_NEAR(vertex_1_on_physical[0], surrogate_sub_model_part_outer.GetNode(iga_model_part.NumberOfGeometries()-current_id).X() , tolerance);
        KRATOS_EXPECT_NEAR(vertex_1_on_physical[1], surrogate_sub_model_part_outer.GetNode(iga_model_part.NumberOfGeometries()-current_id).Y() , tolerance);
    }

    // Ensure the number of nodes matches expectation
    KRATOS_EXPECT_EQ(iga_model_part.Geometries().size(), 10);

}


KRATOS_TEST_CASE_IN_SUITE(TestCreateBrepsSbmUtilitiesInner, KratosIgaFastSuite)
{
    using IndexType = std::size_t;
    using GeometryType = Geometry<Node>;
    using ContainerNodeType = PointerVector<Node>;
    using ContainerEmbeddedNodeType = PointerVector<Point>;
    using BrepCurveOnSurfaceType = BrepCurveOnSurface<ContainerNodeType, true, ContainerEmbeddedNodeType>;
    using CoordinatesArrayType = GeometryType::CoordinatesArrayType;
    
    // Initialization
    std::size_t EchoLevel = 0;
    Model model;
    ModelPart& iga_model_part = model.CreateModelPart("IgaModelPart");
    ModelPart& surrogate_sub_model_part_outer = iga_model_part.CreateSubModelPart("surrogate_outer");
    ModelPart& surrogate_sub_model_part_inner = iga_model_part.CreateSubModelPart("surrogate_inner");
    
    surrogate_sub_model_part_inner.CreateNewProperties(0);
    
    surrogate_sub_model_part_inner.CreateNewNode(1, -1.0, -1.0, 0.0);
    surrogate_sub_model_part_inner.CreateNewNode(2, -1.0, -0.9, 0.0);
    surrogate_sub_model_part_inner.CreateNewNode(3, -1.0, -0.8, 0.0);

    surrogate_sub_model_part_inner.CreateNewNode(4, -0.9, -0.8, 0.0);
    surrogate_sub_model_part_inner.CreateNewNode(5, -0.8, -0.8, 0.0);
    
    surrogate_sub_model_part_inner.CreateNewNode(6, -0.8, -0.9, 0.0);
    surrogate_sub_model_part_inner.CreateNewNode(7, -0.8, -1.0, 0.0);

    surrogate_sub_model_part_inner.CreateNewNode(8, -0.9, -1.0, 0.0);
    
    
    Properties::Pointer p_prop = surrogate_sub_model_part_inner.pGetProperties(0);
    surrogate_sub_model_part_inner.CreateNewCondition("LineCondition2D2N", 1, {{1, 2}}, p_prop);
    surrogate_sub_model_part_inner.CreateNewCondition("LineCondition2D2N", 2, {{2, 3}}, p_prop);
    surrogate_sub_model_part_inner.CreateNewCondition("LineCondition2D2N", 3, {{3, 4}}, p_prop);
    surrogate_sub_model_part_inner.CreateNewCondition("LineCondition2D2N", 4, {{4, 5}}, p_prop);
    surrogate_sub_model_part_inner.CreateNewCondition("LineCondition2D2N", 5, {{5, 6}}, p_prop);
    surrogate_sub_model_part_inner.CreateNewCondition("LineCondition2D2N", 6, {{6, 7}}, p_prop);
    surrogate_sub_model_part_inner.CreateNewCondition("LineCondition2D2N", 7, {{7, 8}}, p_prop);


    // Create "fictituos element" to memorize starting and ending node id for each surrogate boundary loop
    std::vector<ModelPart::IndexType> elem_nodes{1, 8};
    surrogate_sub_model_part_inner.CreateNewElement("Element2D2N", 1, elem_nodes, p_prop);

    Kratos::Parameters parameters(R"(
        {
            "echo_level": 0,
            "model_part_name": "IgaModelPart",
            "geometry_name"  : "nurbs_surface",
            "lower_point_xyz": [-2, -2, 0],
            "upper_point_xyz": [0, 0, 0],
            "lower_point_uvw": [-2, -2, 0],
            "upper_point_uvw": [0, 0, 0],
            "polynomial_order":     [2, 2],
            "number_of_knot_spans": [20, 20]
        }
    )");

    const Point& A_uvw = Point(-2, -2, 0);
    const Point& B_uvw = Point(0, 0, 0);

    NurbsGeometryModeler nurbs_modeler(model, parameters);
    nurbs_modeler.SetupGeometryModel();
    
    NurbsSurfaceGeometryPointerType surface = std::dynamic_pointer_cast<NurbsSurfaceGeometryType>(iga_model_part.pGetGeometry("nurbs_surface"));
    
    // Create the breps for the outer sbm boundary
    CreateBrepsSbmUtilities<Node, Point> CreateBrepsSbmUtilities(EchoLevel);
    CreateBrepsSbmUtilities.CreateSurrogateBoundary(surface, iga_model_part, surrogate_sub_model_part_inner, surrogate_sub_model_part_outer, A_uvw, B_uvw);
    
    const double tolerance = 1e-12;

    for (IndexType current_id = 6; current_id < iga_model_part.NumberOfGeometries(); current_id++) {

        BrepCurveOnSurfaceType::Pointer brep_curve = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(iga_model_part.pGetGeometry(current_id));
        auto ppp = brep_curve->DomainInterval();
        CoordinatesArrayType vertex_1_on_physical;
        CoordinatesArrayType vertex_2_on_physical;
        CoordinatesArrayType local_coord_vertex_1 = ZeroVector(3); 
        CoordinatesArrayType local_coord_vertex_2 = ZeroVector(3);

        local_coord_vertex_1[0] = ppp.GetT0(); local_coord_vertex_2[0] = ppp.GetT1();

        brep_curve->GlobalCoordinates(vertex_1_on_physical, local_coord_vertex_1);
        brep_curve->GlobalCoordinates(vertex_2_on_physical, local_coord_vertex_2);

        KRATOS_EXPECT_NEAR(vertex_1_on_physical[0], surrogate_sub_model_part_inner.GetNode(current_id-5).X() , tolerance);
        KRATOS_EXPECT_NEAR(vertex_1_on_physical[1], surrogate_sub_model_part_inner.GetNode(current_id-5).Y() , tolerance);
    }

    // Ensure the number of nodes matches expectation
    KRATOS_EXPECT_EQ(iga_model_part.Geometries().size(), 14);

}

}