//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    OpenAI
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "custom_processes/snake_sbm_process.h"
#include "custom_modelers/nurbs_geometry_modeler_sbm.h"
#include "custom_modelers/iga_modeler_sbm.h"
#include "iga_application_variables.h"

namespace Kratos::Testing
{
namespace
{
using SizeType = std::size_t;

void CreateCubeOuterSkin(
    ModelPart& rSkinModelPart,
    const double MinCoord = 0.25,
    const double MaxCoord = 1.75)
{
    rSkinModelPart.CreateNewProperties(0);
    auto p_prop = rSkinModelPart.pGetProperties(0);

    rSkinModelPart.CreateNewNode(1, MinCoord, MinCoord, MinCoord);
    rSkinModelPart.CreateNewNode(2, MaxCoord, MinCoord, MinCoord);
    rSkinModelPart.CreateNewNode(3, MaxCoord, MaxCoord, MinCoord);
    rSkinModelPart.CreateNewNode(4, MinCoord, MaxCoord, MinCoord);
    rSkinModelPart.CreateNewNode(5, MinCoord, MinCoord, MaxCoord);
    rSkinModelPart.CreateNewNode(6, MaxCoord, MinCoord, MaxCoord);
    rSkinModelPart.CreateNewNode(7, MaxCoord, MaxCoord, MaxCoord);
    rSkinModelPart.CreateNewNode(8, MinCoord, MaxCoord, MaxCoord);

    // Bottom z = min, outward normal -z
    rSkinModelPart.CreateNewCondition("SurfaceCondition3D3N", 1,  {{1, 3, 2}}, p_prop);
    rSkinModelPart.CreateNewCondition("SurfaceCondition3D3N", 2,  {{1, 4, 3}}, p_prop);
    // Top z = max, outward normal +z
    rSkinModelPart.CreateNewCondition("SurfaceCondition3D3N", 3,  {{5, 6, 7}}, p_prop);
    rSkinModelPart.CreateNewCondition("SurfaceCondition3D3N", 4,  {{5, 7, 8}}, p_prop);
    // Front y = min, outward normal -y
    rSkinModelPart.CreateNewCondition("SurfaceCondition3D3N", 5,  {{1, 2, 6}}, p_prop);
    rSkinModelPart.CreateNewCondition("SurfaceCondition3D3N", 6,  {{1, 6, 5}}, p_prop);
    // Back y = max, outward normal +y
    rSkinModelPart.CreateNewCondition("SurfaceCondition3D3N", 7,  {{4, 8, 7}}, p_prop);
    rSkinModelPart.CreateNewCondition("SurfaceCondition3D3N", 8,  {{4, 7, 3}}, p_prop);
    // Left x = min, outward normal -x
    rSkinModelPart.CreateNewCondition("SurfaceCondition3D3N", 9,  {{1, 5, 8}}, p_prop);
    rSkinModelPart.CreateNewCondition("SurfaceCondition3D3N", 10, {{1, 8, 4}}, p_prop);
    // Right x = max, outward normal +x
    rSkinModelPart.CreateNewCondition("SurfaceCondition3D3N", 11, {{2, 7, 6}}, p_prop);
    rSkinModelPart.CreateNewCondition("SurfaceCondition3D3N", 12, {{2, 3, 7}}, p_prop);
}

void CreateOctahedronInnerSkin(
    ModelPart& rSkinModelPart,
    const double Center = 1.0,
    const double Radius = 0.35)
{
    rSkinModelPart.CreateNewProperties(0);
    auto p_prop = rSkinModelPart.pGetProperties(0);

    rSkinModelPart.CreateNewNode(1, Center + Radius, Center, Center);
    rSkinModelPart.CreateNewNode(2, Center - Radius, Center, Center);
    rSkinModelPart.CreateNewNode(3, Center, Center + Radius, Center);
    rSkinModelPart.CreateNewNode(4, Center, Center - Radius, Center);
    rSkinModelPart.CreateNewNode(5, Center, Center, Center + Radius);
    rSkinModelPart.CreateNewNode(6, Center, Center, Center - Radius);

    rSkinModelPart.CreateNewCondition("SurfaceCondition3D3N", 1, {{1, 3, 5}}, p_prop);
    rSkinModelPart.CreateNewCondition("SurfaceCondition3D3N", 2, {{3, 2, 5}}, p_prop);
    rSkinModelPart.CreateNewCondition("SurfaceCondition3D3N", 3, {{2, 4, 5}}, p_prop);
    rSkinModelPart.CreateNewCondition("SurfaceCondition3D3N", 4, {{4, 1, 5}}, p_prop);
    rSkinModelPart.CreateNewCondition("SurfaceCondition3D3N", 5, {{3, 1, 6}}, p_prop);
    rSkinModelPart.CreateNewCondition("SurfaceCondition3D3N", 6, {{2, 3, 6}}, p_prop);
    rSkinModelPart.CreateNewCondition("SurfaceCondition3D3N", 7, {{4, 2, 6}}, p_prop);
    rSkinModelPart.CreateNewCondition("SurfaceCondition3D3N", 8, {{1, 4, 6}}, p_prop);
}

Vector MakeVector(const std::vector<double>& rValues)
{
    Vector values(rValues.size());
    for (IndexType i = 0; i < rValues.size(); ++i) {
        values[i] = rValues[i];
    }
    return values;
}

bool IsOnBoundaryPlane(const Point& rCenter, const double Tolerance)
{
    return std::abs(rCenter.X()) < Tolerance
        || std::abs(rCenter.X() - 2.0) < Tolerance
        || std::abs(rCenter.Y()) < Tolerance
        || std::abs(rCenter.Y() - 2.0) < Tolerance
        || std::abs(rCenter.Z()) < Tolerance
        || std::abs(rCenter.Z() - 2.0) < Tolerance;
}

double ExpectedOuterOrientationSign(const Point& rCenter, const double Tolerance)
{
    const bool is_on_box_boundary =
        std::abs(rCenter.X()) < Tolerance || std::abs(rCenter.X() - 2.0) < Tolerance ||
        std::abs(rCenter.Y()) < Tolerance || std::abs(rCenter.Y() - 2.0) < Tolerance ||
        std::abs(rCenter.Z()) < Tolerance || std::abs(rCenter.Z() - 2.0) < Tolerance;

    if (is_on_box_boundary) {
        return 1.0;
    }

    const bool is_on_inner_cube_boundary =
        std::abs(rCenter.X() - 0.5) < Tolerance || std::abs(rCenter.X() - 1.5) < Tolerance ||
        std::abs(rCenter.Y() - 0.5) < Tolerance || std::abs(rCenter.Y() - 1.5) < Tolerance ||
        std::abs(rCenter.Z() - 0.5) < Tolerance || std::abs(rCenter.Z() - 1.5) < Tolerance;

    if (is_on_inner_cube_boundary) {
        return -1.0;
    }

    KRATOS_ERROR << "Support outer condition center is not on an expected plane: " << rCenter << std::endl;
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(SnakeSbmProcessCubeOuter3D, KratosIgaFastSuite)
{
    Model model;
    auto& r_iga_model_part = model.CreateModelPart("iga_model_part");
    r_iga_model_part.CreateSubModelPart("surrogate_inner");
    auto& r_surrogate_outer = r_iga_model_part.CreateSubModelPart("surrogate_outer");

    model.CreateModelPart("skin_model_part_inner_initial");
    auto& r_skin_outer_initial = model.CreateModelPart("skin_model_part_outer_initial");
    CreateCubeOuterSkin(r_skin_outer_initial);

    auto& r_skin_model_part = model.CreateModelPart("skin_model_part");
    r_skin_model_part.CreateSubModelPart("inner");
    r_skin_model_part.CreateSubModelPart("outer");

    r_iga_model_part.SetValue(KNOT_VECTOR_U, MakeVector({0.0, 1.0, 2.0}));
    r_iga_model_part.SetValue(KNOT_VECTOR_V, MakeVector({0.0, 1.0, 2.0}));
    r_iga_model_part.SetValue(KNOT_VECTOR_W, MakeVector({0.0, 1.0, 2.0}));

    Parameters snake_parameters(R"(
        {
            "model_part_name" : "iga_model_part",
            "skin_model_part_inner_initial_name" : "skin_model_part_inner_initial",
            "skin_model_part_outer_initial_name" : "skin_model_part_outer_initial",
            "skin_model_part_name" : "skin_model_part",
            "echo_level" : 0,
            "lambda_inner" : 0.5,
            "lambda_outer" : 0.5,
            "number_of_inner_loops": 0
        }
    )");

    SnakeSbmProcess snake_sbm_process(model, snake_parameters);
    snake_sbm_process.Execute();

    const double tolerance = 1.0e-10;
    KRATOS_EXPECT_EQ(r_surrogate_outer.NumberOfConditions(), 12);
    KRATOS_EXPECT_EQ(r_surrogate_outer.NumberOfNodes(), 27);

    for (const auto& r_condition : r_surrogate_outer.Conditions()) {
        KRATOS_EXPECT_EQ(r_condition.GetGeometry().PointsNumber(), 4);
        KRATOS_EXPECT_TRUE(IsOnBoundaryPlane(r_condition.GetGeometry().Center(), tolerance));
    }
}

KRATOS_TEST_CASE_IN_SUITE(SnakeSbmProcessOctahedronInner3D, KratosIgaFastSuite)
{
    Model model;
    auto& r_iga_model_part = model.CreateModelPart("iga_model_part");
    auto& r_surrogate_inner = r_iga_model_part.CreateSubModelPart("surrogate_inner");
    r_iga_model_part.CreateSubModelPart("surrogate_outer");

    auto& r_skin_inner_initial = model.CreateModelPart("skin_model_part_inner_initial");
    CreateOctahedronInnerSkin(r_skin_inner_initial);
    model.CreateModelPart("skin_model_part_outer_initial");

    auto& r_skin_model_part = model.CreateModelPart("skin_model_part");
    r_skin_model_part.CreateSubModelPart("inner");
    r_skin_model_part.CreateSubModelPart("outer");

    r_iga_model_part.SetValue(KNOT_VECTOR_U, MakeVector({0.0, 0.5, 1.0, 1.5, 2.0}));
    r_iga_model_part.SetValue(KNOT_VECTOR_V, MakeVector({0.0, 0.5, 1.0, 1.5, 2.0}));
    r_iga_model_part.SetValue(KNOT_VECTOR_W, MakeVector({0.0, 0.5, 1.0, 1.5, 2.0}));

    Parameters snake_parameters(R"(
        {
            "model_part_name" : "iga_model_part",
            "skin_model_part_inner_initial_name" : "skin_model_part_inner_initial",
            "skin_model_part_outer_initial_name" : "skin_model_part_outer_initial",
            "skin_model_part_name" : "skin_model_part",
            "echo_level" : 0,
            "lambda_inner" : 0.5,
            "lambda_outer" : 0.5,
            "number_of_inner_loops": 1
        }
    )");

    SnakeSbmProcess snake_sbm_process(model, snake_parameters);
    snake_sbm_process.Execute();

    const double tolerance = 1.0e-10;
    KRATOS_EXPECT_EQ(r_surrogate_inner.NumberOfConditions(), 48);
    KRATOS_EXPECT_EQ(r_surrogate_inner.NumberOfNodes(), 125);
    KRATOS_EXPECT_EQ(r_surrogate_inner.NumberOfElements(), 1);

    IndexType boundary_true_count = 0;
    IndexType boundary_false_count = 0;
    for (const auto& r_condition : r_surrogate_inner.Conditions()) {
        KRATOS_EXPECT_EQ(r_condition.GetGeometry().PointsNumber(), 4);
        const auto center = r_condition.GetGeometry().Center();
        KRATOS_EXPECT_GE(center.X(), -tolerance);
        KRATOS_EXPECT_LE(center.X(), 2.0 + tolerance);
        KRATOS_EXPECT_GE(center.Y(), -tolerance);
        KRATOS_EXPECT_LE(center.Y(), 2.0 + tolerance);
        KRATOS_EXPECT_GE(center.Z(), -tolerance);
        KRATOS_EXPECT_LE(center.Z(), 2.0 + tolerance);
        if (r_condition.Is(BOUNDARY)) {
            ++boundary_true_count;
        } else {
            ++boundary_false_count;
        }
    }

    KRATOS_EXPECT_EQ(boundary_true_count, 48);
    KRATOS_EXPECT_EQ(boundary_false_count, 0);

    const auto& r_loop_geometry = r_surrogate_inner.ElementsBegin()->GetGeometry();
    KRATOS_EXPECT_EQ(r_loop_geometry.PointsNumber(), 2);
    KRATOS_EXPECT_EQ(r_loop_geometry[0].Id(), 1);
    KRATOS_EXPECT_EQ(r_loop_geometry[1].Id(), 48);
}

KRATOS_TEST_CASE_IN_SUITE(BrepVolumeQuadraturePointGenerationOuter3D, KratosIgaFastSuite)
{
    Model model;
    auto& r_iga_model_part = model.CreateModelPart("IgaModelPart");
    r_iga_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

    auto& r_skin_outer_initial = model.CreateModelPart("skin_model_part_outer_initial");
    CreateCubeOuterSkin(r_skin_outer_initial, 0.5, 1.5);

    Parameters nurbs_modeler_parameters(R"(
        {
            "model_part_name" : "IgaModelPart",
            "lower_point_xyz": [0.0, 0.0, 0.0],
            "upper_point_xyz": [2.0, 2.0, 2.0],
            "lower_point_uvw": [0.0, 0.0, 0.0],
            "upper_point_uvw": [2.0, 2.0, 2.0],
            "polynomial_order" : [1, 1, 1],
            "number_of_knot_spans" : [4, 4, 4],
            "lambda_outer": 0.5,
            "number_of_inner_loops": 0,
            "skin_model_part_outer_initial_name": "skin_model_part_outer_initial",
            "skin_model_part_name": "skin_model_part",
            "echo_level": 0
        }
    )");

    NurbsGeometryModelerSbm nurbs_modeler(model, nurbs_modeler_parameters);
    nurbs_modeler.SetupGeometryModel();
    nurbs_modeler.PrepareGeometryModel();
    nurbs_modeler.SetupModelPart();

    auto& r_brep_volume = r_iga_model_part.GetGeometry(1);
    Geometry<Node>::GeometriesArrayType quadrature_geometries;
    auto integration_info = r_brep_volume.GetDefaultIntegrationInfo();
    r_brep_volume.CreateQuadraturePointGeometries(quadrature_geometries, 2, integration_info);

    KRATOS_EXPECT_EQ(r_iga_model_part.NumberOfGeometries(), 73);
    KRATOS_EXPECT_EQ(model.GetModelPart("IgaModelPart.surrogate_outer").NumberOfConditions(), 72);
    KRATOS_EXPECT_EQ(quadrature_geometries.size(), 64);
}

KRATOS_TEST_CASE_IN_SUITE(IgaModelerSbmSupportOuter3D, KratosIgaFastSuite)
{
    Model model;
    auto& r_iga_model_part = model.CreateModelPart("IgaModelPart");
    r_iga_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

    auto& r_skin_outer_initial = model.CreateModelPart("skin_model_part_outer_initial");
    CreateCubeOuterSkin(r_skin_outer_initial, 0.5, 1.5);

    Parameters nurbs_modeler_parameters(R"(
        {
            "model_part_name" : "IgaModelPart",
            "lower_point_xyz": [0.0, 0.0, 0.0],
            "upper_point_xyz": [2.0, 2.0, 2.0],
            "lower_point_uvw": [0.0, 0.0, 0.0],
            "upper_point_uvw": [2.0, 2.0, 2.0],
            "polynomial_order" : [1, 1, 1],
            "number_of_knot_spans" : [4, 4, 4],
            "lambda_outer": 0.5,
            "number_of_inner_loops": 0,
            "skin_model_part_outer_initial_name": "skin_model_part_outer_initial",
            "skin_model_part_name": "skin_model_part",
            "echo_level": 0
        }
    )");

    NurbsGeometryModelerSbm nurbs_modeler(model, nurbs_modeler_parameters);
    nurbs_modeler.SetupGeometryModel();
    nurbs_modeler.PrepareGeometryModel();
    nurbs_modeler.SetupModelPart();

    Parameters iga_modeler_parameters(R"(
        {
            "echo_level": 0,
            "skin_model_part_name": "skin_model_part",
            "analysis_model_part_name": "IgaModelPart",
            "element_condition_list": [
                {
                    "geometry_type": "SurfaceEdge",
                    "iga_model_part": "SBM_Support_outer",
                    "type": "condition",
                    "name": "SbmLaplacianConditionDirichlet",
                    "shape_function_derivatives_order": 2,
                    "sbm_parameters": {
                        "is_inner" : false
                    }
                }
            ]
        }
    )");

    IgaModelerSbm iga_modeler(model, iga_modeler_parameters);
    iga_modeler.SetupGeometryModel();
    iga_modeler.PrepareGeometryModel();
    iga_modeler.SetupModelPart();

    auto& r_support_model_part = model.GetModelPart("IgaModelPart.SBM_Support_outer");

    KRATOS_EXPECT_FALSE(model.HasModelPart("IgaModelPart.ConvectionDiffusionDomain"));
    KRATOS_EXPECT_EQ(r_support_model_part.NumberOfConditions(), 288);
    KRATOS_EXPECT_EQ(r_support_model_part.NumberOfNodes(), 0);

    const auto& r_condition = *(r_support_model_part.ConditionsBegin());
    KRATOS_EXPECT_EQ(r_condition.GetValue(IDENTIFIER), "outer");
    KRATOS_EXPECT_EQ(r_condition.GetValue(NEIGHBOUR_CONDITIONS).size(), 1);

    const auto& r_condition_knot_span_sizes = r_condition.GetValue(KNOT_SPAN_SIZES);
    KRATOS_EXPECT_EQ(r_condition_knot_span_sizes.size(), 3);
    KRATOS_EXPECT_NEAR(r_condition_knot_span_sizes[0], 0.5, 1.0e-12);
    KRATOS_EXPECT_NEAR(r_condition_knot_span_sizes[1], 0.5, 1.0e-12);
    KRATOS_EXPECT_NEAR(r_condition_knot_span_sizes[2], 0.5, 1.0e-12);

    const double normal_tolerance = 1.0e-10;
    array_1d<double, 3> domain_center;
    domain_center[0] = 1.0;
    domain_center[1] = 1.0;
    domain_center[2] = 1.0;
    for (const auto& r_support_condition : r_support_model_part.Conditions()) {
        array_1d<double, 3> computed_normal;
        r_support_condition.GetGeometry().Calculate(NORMAL, computed_normal);
        computed_normal /= norm_2(computed_normal);
        KRATOS_EXPECT_NEAR(norm_2(computed_normal), 1.0, normal_tolerance);

        const array_1d<double, 3> radial_direction = r_support_condition.GetGeometry().Center() - domain_center;
        KRATOS_EXPECT_GT(std::abs(inner_prod(computed_normal, radial_direction)), 0.0);
        KRATOS_EXPECT_NEAR(
            inner_prod(computed_normal, radial_direction),
            ExpectedOuterOrientationSign(r_support_condition.GetGeometry().Center(), normal_tolerance) * std::abs(inner_prod(computed_normal, radial_direction)),
            normal_tolerance);
    }
}

KRATOS_TEST_CASE_IN_SUITE(IgaModelerSbmSupportInner3D, KratosIgaFastSuite)
{
    Model model;
    auto& r_iga_model_part = model.CreateModelPart("IgaModelPart");
    r_iga_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

    auto& r_skin_inner_initial = model.CreateModelPart("skin_model_part_inner_initial");
    CreateOctahedronInnerSkin(r_skin_inner_initial);
    model.CreateModelPart("skin_model_part_outer_initial");

    Parameters nurbs_modeler_parameters(R"(
        {
            "model_part_name" : "IgaModelPart",
            "lower_point_xyz": [0.0, 0.0, 0.0],
            "upper_point_xyz": [2.0, 2.0, 2.0],
            "lower_point_uvw": [0.0, 0.0, 0.0],
            "upper_point_uvw": [2.0, 2.0, 2.0],
            "polynomial_order" : [1, 1, 1],
            "number_of_knot_spans" : [4, 4, 4],
            "lambda_inner": 0.5,
            "number_of_inner_loops": 1,
            "skin_model_part_inner_initial_name": "skin_model_part_inner_initial",
            "skin_model_part_name": "skin_model_part",
            "echo_level": 0
        }
    )");

    NurbsGeometryModelerSbm nurbs_modeler(model, nurbs_modeler_parameters);
    nurbs_modeler.SetupGeometryModel();
    nurbs_modeler.PrepareGeometryModel();
    nurbs_modeler.SetupModelPart();

    Parameters iga_modeler_parameters(R"(
        {
            "echo_level": 0,
            "skin_model_part_name": "skin_model_part",
            "analysis_model_part_name": "IgaModelPart",
            "element_condition_list": [
                {
                    "geometry_type": "SurfaceEdge",
                    "iga_model_part": "SBM_Support_inner",
                    "type": "condition",
                    "name": "SbmLaplacianConditionDirichlet",
                    "shape_function_derivatives_order": 2,
                    "sbm_parameters": {
                        "is_inner" : true
                    }
                }
            ]
        }
    )");

    IgaModelerSbm iga_modeler(model, iga_modeler_parameters);
    iga_modeler.SetupGeometryModel();
    iga_modeler.PrepareGeometryModel();
    iga_modeler.SetupModelPart();

    auto& r_surrogate_inner = model.GetModelPart("IgaModelPart.surrogate_inner");
    auto& r_support_model_part = model.GetModelPart("IgaModelPart.SBM_Support_inner");

    KRATOS_EXPECT_EQ(r_surrogate_inner.NumberOfConditions(), 48);
    KRATOS_EXPECT_EQ(r_surrogate_inner.NumberOfNodes(), 125);
    KRATOS_EXPECT_EQ(r_surrogate_inner.NumberOfElements(), 1);
    KRATOS_EXPECT_EQ(r_support_model_part.NumberOfConditions(), 192);
    KRATOS_EXPECT_EQ(r_support_model_part.NumberOfNodes(), 0);

    const auto& r_condition = *(r_support_model_part.ConditionsBegin());
    KRATOS_EXPECT_EQ(r_condition.GetValue(IDENTIFIER), "inner");
    KRATOS_EXPECT_EQ(r_condition.GetValue(NEIGHBOUR_CONDITIONS).size(), 1);

    const auto& r_condition_knot_span_sizes = r_condition.GetValue(KNOT_SPAN_SIZES);
    KRATOS_EXPECT_EQ(r_condition_knot_span_sizes.size(), 3);
    KRATOS_EXPECT_NEAR(r_condition_knot_span_sizes[0], 0.5, 1.0e-12);
    KRATOS_EXPECT_NEAR(r_condition_knot_span_sizes[1], 0.5, 1.0e-12);
    KRATOS_EXPECT_NEAR(r_condition_knot_span_sizes[2], 0.5, 1.0e-12);

    const double normal_tolerance = 1.0e-10;
    array_1d<double, 3> domain_center;
    domain_center[0] = 1.0;
    domain_center[1] = 1.0;
    domain_center[2] = 1.0;
    for (const auto& r_support_condition : r_support_model_part.Conditions()) {
        array_1d<double, 3> computed_normal;
        r_support_condition.GetGeometry().Calculate(NORMAL, computed_normal);
        computed_normal /= norm_2(computed_normal);
        KRATOS_EXPECT_NEAR(norm_2(computed_normal), 1.0, normal_tolerance);

        const array_1d<double, 3> radial_direction = r_support_condition.GetGeometry().Center() - domain_center;
        KRATOS_EXPECT_LT(inner_prod(computed_normal, radial_direction), 0.0);
    }
}
} // namespace Kratos::Testing
