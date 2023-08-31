//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// Project includes
#include "testing/testing.h"
#include "geometries/oriented_bounding_box.h"
#include "includes/expect.h"
#include "geometries/point.h"
#include "utilities/intersection_utilities.h"

// // Debug
// #include "containers/model.h"
// #include "includes/gid_io.h"

namespace Kratos {
    namespace Testing {

    /**
     *  Here we check the HasIntersection function of the OBB 2D
     * We compare with the intersection of two quadrilaterals
     */
    KRATOS_TEST_CASE_IN_SUITE(OBBHasIntersection2D, KratosCoreFastSuite)
    {
        array_1d<double, 3> first_center;
        first_center[0] = 0.5;
        first_center[1] = 0.5;
        first_center[2] = 0.0;

        array_1d<array_1d<double, 3>, 2> first_directions;
        first_directions[0][0] = std::sqrt(2.0)/2.0;
        first_directions[0][1] = std::sqrt(2.0)/2.0;
        first_directions[0][2] = 0.0;
        first_directions[1][0] = -std::sqrt(2.0)/2.0;
        first_directions[1][1] = std::sqrt(2.0)/2.0;
        first_directions[1][2] = 0.0;

        array_1d<double, 2> first_half_lenghts;
        first_half_lenghts[0] = std::sqrt(2.0)/2.0;
        first_half_lenghts[1] = 0.1;

        OrientedBoundingBox<2> first_obb(first_center, first_directions, first_half_lenghts);

        array_1d<double, 3> second_center;
        second_center[0] = 0.0;
        second_center[1] = 0.0;
        second_center[2] = 0.0;

        array_1d<array_1d<double, 3>, 2> second_directions;
        second_directions[0][0] = 1.0;
        second_directions[0][1] = 0.0;
        second_directions[0][2] = 0.0;
        second_directions[1][0] = 0.0;
        second_directions[1][1] = 1.0;
        second_directions[1][2] = 0.0;

        array_1d<double, 2> second_half_lenghts;
        second_half_lenghts[0] = 0.5;
        second_half_lenghts[1] = 0.1;

        OrientedBoundingBox<2> second_obb(second_center, second_directions, second_half_lenghts);

        auto first_quad = first_obb.GetEquivalentGeometry();
        auto second_quad = second_obb.GetEquivalentGeometry();

        bool has_intersection_reference = false;

        array_1d<double, 3> local_coords;
        for (auto& r_point : first_quad) {
            if (second_quad.IsInside(r_point.Coordinates(), local_coords)) {
                has_intersection_reference = true;
                break;
            }
        }
        for (auto& r_point : second_quad) {
            if (first_quad.IsInside(r_point.Coordinates(), local_coords)) {
                has_intersection_reference = true;
                break;
            }
        }

//         // Debug
//         first_obb.GetEquivalentRotatedGeometry(second_quad);
//         first_obb.GetEquivalentRotatedGeometry(first_quad);
//
//         Model current_model;
//         ModelPart& r_main_model_part = current_model.CreateModelPart("Main");
//         ModelPart& r_model_part1 = r_main_model_part.CreateSubModelPart("1");
//         ModelPart& r_model_part2 = r_main_model_part.CreateSubModelPart("2");
//
//         for (int i = 0; i < 4; ++i) {
//             r_model_part1.CreateNewNode(i + 1, first_quad[i].X(), first_quad[i].Y(), first_quad[i].Z());
//             r_model_part2.CreateNewNode(i + 5, second_quad[i].X(), second_quad[i].Y(), second_quad[i].Z());
//         }
//
//         Properties::Pointer p_properties_0 = Kratos::make_shared<Properties>(0);
//         Properties::Pointer p_properties_1 = Kratos::make_shared<Properties>(1);
//
//         r_model_part1.CreateNewElement("Element2D4N", 1, {1, 2, 3, 4}, p_properties_0);
//         r_model_part2.CreateNewElement("Element2D4N", 2, {5, 6, 7, 8}, p_properties_1);
//
//         GidIO<> gid_io("test", GiD_PostBinary, SingleFile, WriteDeformed, WriteConditions);
//         gid_io.InitializeMesh(0.0);
//         gid_io.WriteMesh(r_main_model_part.GetMesh());
//         gid_io.FinalizeMesh();
//         gid_io.InitializeResults(0, r_main_model_part.GetMesh());
//         gid_io.FinalizeResults();

        bool has_intersection = first_obb.HasIntersection(second_obb);

        KRATOS_EXPECT_EQ(has_intersection_reference, has_intersection);

        // Moving up tyhe second OBBHasIntersection2D
        second_center[1] = 0.5;
        second_obb.SetCenter(second_center);
        second_quad = second_obb.GetEquivalentGeometry();

//         // Debug
//         first_obb.GetEquivalentRotatedGeometry(second_quad);
//         first_obb.GetEquivalentRotatedGeometry(first_quad);
//
//         Model current_model;
//         ModelPart& r_main_model_part = current_model.CreateModelPart("Main");
//         ModelPart& r_model_part1 = r_main_model_part.CreateSubModelPart("1");
//         ModelPart& r_model_part2 = r_main_model_part.CreateSubModelPart("2");
//
//         for (int i = 0; i < 4; ++i) {
//             r_model_part1.CreateNewNode(i + 1, first_quad[i].X(), first_quad[i].Y(), first_quad[i].Z());
//             r_model_part2.CreateNewNode(i + 5, second_quad[i].X(), second_quad[i].Y(), second_quad[i].Z());
//         }
//
//         Properties::Pointer p_properties_0 = Kratos::make_shared<Properties>(0);
//         Properties::Pointer p_properties_1 = Kratos::make_shared<Properties>(1);
//
//         r_model_part1.CreateNewElement("Element2D4N", 1, {1, 2, 3, 4}, p_properties_0);
//         r_model_part2.CreateNewElement("Element2D4N", 2, {5, 6, 7, 8}, p_properties_1);
//
//         GidIO<> gid_io("test", GiD_PostBinary, SingleFile, WriteDeformed, WriteConditions);
//         gid_io.InitializeMesh(0.0);
//         gid_io.WriteMesh(r_main_model_part.GetMesh());
//         gid_io.FinalizeMesh();
//         gid_io.InitializeResults(0, r_main_model_part.GetMesh());
//         gid_io.FinalizeResults();

        has_intersection_reference = false;

        for (auto& r_point : first_quad) {
            if (second_quad.IsInside(r_point.Coordinates(), local_coords)) {
                has_intersection_reference = true;
                break;
            }
        }
        for (auto& r_point : second_quad) {
            if (first_quad.IsInside(r_point.Coordinates(), local_coords)) {
                has_intersection_reference = true;
                break;
            }
        }

        has_intersection = first_obb.HasIntersection(second_obb);

        KRATOS_EXPECT_EQ(has_intersection_reference, has_intersection);
    }

    /**
     *  Here we check the HasIntersection function of the OBB 2D
     * We compare with the intersection of two quadrilaterals
     */
    KRATOS_TEST_CASE_IN_SUITE(OBBHasIntersection2DCrossed, KratosCoreFastSuite)
    {
        array_1d<double, 3> first_center;
        first_center[0] = 0.0;
        first_center[1] = 0.0;
        first_center[2] = 0.0;

        array_1d<array_1d<double, 3>, 2> first_directions;
        first_directions[0][0] = 1.0;
        first_directions[0][1] = 0.0;
        first_directions[0][2] = 0.0;
        first_directions[1][0] = 0.0;
        first_directions[1][1] = 1.0;
        first_directions[1][2] = 0.0;

        array_1d<double, 2> first_half_lenghts;
        first_half_lenghts[0] = 0.5;
        first_half_lenghts[1] = 0.1;

        OrientedBoundingBox<2> first_obb(first_center, first_directions, first_half_lenghts);

        array_1d<double, 3> second_center;
        second_center[0] = 0.0;
        second_center[1] = 0.0;
        second_center[2] = 0.0;

        array_1d<array_1d<double, 3>, 2> second_directions;
        second_directions[0][0] = 0.0;
        second_directions[0][1] = 1.0;
        second_directions[0][2] = 0.0;
        second_directions[1][0] = 1.0;
        second_directions[1][1] = 0.0;
        second_directions[1][2] = 0.0;

        array_1d<double, 2> second_half_lenghts;
        second_half_lenghts[0] = 0.5;
        second_half_lenghts[1] = 0.1;

        OrientedBoundingBox<2> second_obb(second_center, second_directions, second_half_lenghts);

        auto first_quad = first_obb.GetEquivalentGeometry();
        auto second_quad = second_obb.GetEquivalentGeometry();

        bool has_intersection_reference = false;

        array_1d<double, 3> local_coords;
        for (auto& r_point : first_quad) {
            if (second_quad.IsInside(r_point.Coordinates(), local_coords)) {
                has_intersection_reference = true;
                break;
            }
        }
        for (auto& r_point : second_quad) {
            if (first_quad.IsInside(r_point.Coordinates(), local_coords)) {
                has_intersection_reference = true;
                break;
            }
        }

        auto r_edges_1 = first_quad.GenerateEdges();
        auto r_edges_2 = second_quad.GenerateEdges();
        Point int_pt(0.0,0.0,0.0);
        for (auto& r_edge_1 : r_edges_1) {
            for (auto& r_edge_2 : r_edges_2) {
                const int int_id = IntersectionUtilities::ComputeLineLineIntersection(
                    r_edge_1,
                    r_edge_2[0].Coordinates(),
                    r_edge_2[1].Coordinates(),
                    int_pt.Coordinates());

                if (int_id != 0){
                    has_intersection_reference = true;
                    break;
                }
            }
        }

//         // Debug
//         first_obb.GetEquivalentRotatedGeometry(second_quad);
//         first_obb.GetEquivalentRotatedGeometry(first_quad);
//
//         Model current_model;
//         ModelPart& r_main_model_part = current_model.CreateModelPart("Main");
//         ModelPart& r_model_part1 = r_main_model_part.CreateSubModelPart("1");
//         ModelPart& r_model_part2 = r_main_model_part.CreateSubModelPart("2");
//
//         for (int i = 0; i < 4; ++i) {
//             r_model_part1.CreateNewNode(i + 1, first_quad[i].X(), first_quad[i].Y(), first_quad[i].Z());
//             r_model_part2.CreateNewNode(i + 5, second_quad[i].X(), second_quad[i].Y(), second_quad[i].Z());
//         }
//
//         Properties::Pointer p_properties_0 = Kratos::make_shared<Properties>(0);
//         Properties::Pointer p_properties_1 = Kratos::make_shared<Properties>(1);
//
//         r_model_part1.CreateNewElement("Element2D4N", 1, {1, 2, 3, 4}, p_properties_0);
//         r_model_part2.CreateNewElement("Element2D4N", 2, {5, 6, 7, 8}, p_properties_1);
//
//         GidIO<> gid_io("test", GiD_PostBinary, SingleFile, WriteDeformed, WriteConditions);
//         gid_io.InitializeMesh(0.0);
//         gid_io.WriteMesh(r_main_model_part.GetMesh());
//         gid_io.FinalizeMesh();
//         gid_io.InitializeResults(0, r_main_model_part.GetMesh());
//         gid_io.FinalizeResults();

        bool has_intersection = first_obb.HasIntersection(second_obb);

        KRATOS_EXPECT_EQ(has_intersection_reference, has_intersection);
    }

    /**
     *  Here we check the HasIntersection function of the OBB 3D
     * We compare with the intersection of two quadrilaterals
     */
    KRATOS_TEST_CASE_IN_SUITE(OBBHasIntersection3DSimple, KratosCoreFastSuite)
    {
        array_1d<double, 3> first_center;
        first_center[0] = 0.5;
        first_center[1] = 0.5;
        first_center[2] = 0.0;

        array_1d<array_1d<double, 3>, 3> first_directions;
        first_directions[0][0] = std::sqrt(2.0)/2.0;
        first_directions[0][1] = std::sqrt(2.0)/2.0;
        first_directions[0][2] = 0.0;
        first_directions[1][0] = -std::sqrt(2.0)/2.0;
        first_directions[1][1] = std::sqrt(2.0)/2.0;
        first_directions[1][2] = 0.0;
        first_directions[2][0] = 0.0;
        first_directions[2][1] = 0.0;
        first_directions[2][2] = 1.0;

        array_1d<double, 3> first_half_lenghts;
        first_half_lenghts[0] = std::sqrt(2.0)/2.0;
        first_half_lenghts[1] = 0.1;
        first_half_lenghts[2] = 0.1;

        OrientedBoundingBox<3> first_obb(first_center, first_directions, first_half_lenghts);

        array_1d<double, 3> second_center;
        second_center[0] = 0.0;
        second_center[1] = 0.0;
        second_center[2] = 0.0;

        array_1d<array_1d<double, 3>, 3> second_directions;
        second_directions[0][0] = 1.0;
        second_directions[0][1] = 0.0;
        second_directions[0][2] = 0.0;
        second_directions[1][0] = 0.0;
        second_directions[1][1] = 1.0;
        second_directions[1][2] = 0.0;
        second_directions[2][0] = 0.0;
        second_directions[2][1] = 0.0;
        second_directions[2][2] = 1.0;

        array_1d<double, 3> second_half_lenghts;
        second_half_lenghts[0] = 0.5;
        second_half_lenghts[1] = 0.1;
        second_half_lenghts[2] = 0.1;

        OrientedBoundingBox<3> second_obb(second_center, second_directions, second_half_lenghts);

        auto first_hexa = first_obb.GetEquivalentGeometry();
        auto second_hexa = second_obb.GetEquivalentGeometry();

        bool has_intersection_reference = false;

        array_1d<double, 3> local_coords;
        for (auto& r_point : first_hexa) {
            if (second_hexa.IsInside(r_point.Coordinates(), local_coords)) {
                has_intersection_reference = true;
                break;
            }
        }
        for (auto& r_point : second_hexa) {
            if (first_hexa.IsInside(r_point.Coordinates(), local_coords)) {
                has_intersection_reference = true;
                break;
            }
        }

//         // Debug
//         first_obb.GetEquivalentRotatedGeometry(second_hexa);
//         first_obb.GetEquivalentRotatedGeometry(first_hexa);
//
//         Model current_model;
//         ModelPart& r_main_model_part = current_model.CreateModelPart("Main");
//         ModelPart& r_model_part1 = r_main_model_part.CreateSubModelPart("1");
//         ModelPart& r_model_part2 = r_main_model_part.CreateSubModelPart("2");
//
//         for (int i = 0; i < 8; ++i) {
//             r_model_part1.CreateNewNode(i + 1, first_hexa[i].X(), first_hexa[i].Y(), first_hexa[i].Z());
//             r_model_part2.CreateNewNode(i + 9, second_hexa[i].X(), second_hexa[i].Y(), second_hexa[i].Z());
//         }
//
//         Properties::Pointer p_properties_0 = Kratos::make_shared<Properties>(0);
//         Properties::Pointer p_properties_1 = Kratos::make_shared<Properties>(1);
//
//         r_model_part1.CreateNewElement("Element3D8N", 1, {1, 2, 3, 4, 5, 6, 7, 8}, p_properties_0);
//         r_model_part2.CreateNewElement("Element3D8N", 2, {9, 10, 11, 12, 13, 14, 15, 16}, p_properties_1);
//
//         GidIO<> gid_io("test", GiD_PostBinary, SingleFile, WriteDeformed, WriteConditions);
//         gid_io.InitializeMesh(0.0);
//         gid_io.WriteMesh(r_main_model_part.GetMesh());
//         gid_io.FinalizeMesh();
//         gid_io.InitializeResults(0, r_main_model_part.GetMesh());
//         gid_io.FinalizeResults();

        bool has_intersection = first_obb.HasIntersection(second_obb);

        KRATOS_EXPECT_EQ(has_intersection_reference, has_intersection);

        // Moving up tyhe second OBB
        second_center[1] = 0.5;
        second_obb.SetCenter(second_center);
        second_hexa = second_obb.GetEquivalentGeometry();

//         // Debug
//         first_obb.GetEquivalentRotatedGeometry(second_hexa);
//         first_obb.GetEquivalentRotatedGeometry(first_hexa);
//
//         Model current_model;
//         ModelPart& r_main_model_part = current_model.CreateModelPart("Main");
//         ModelPart& r_model_part1 = r_main_model_part.CreateSubModelPart("1");
//         ModelPart& r_model_part2 = r_main_model_part.CreateSubModelPart("2");
//
//         for (int i = 0; i < 8; ++i) {
//             r_model_part1.CreateNewNode(i + 1, first_hexa[i].X(), first_hexa[i].Y(), first_hexa[i].Z());
//             r_model_part2.CreateNewNode(i + 9, second_hexa[i].X(), second_hexa[i].Y(), second_hexa[i].Z());
//         }
//
//         Properties::Pointer p_properties_0 = Kratos::make_shared<Properties>(0);
//         Properties::Pointer p_properties_1 = Kratos::make_shared<Properties>(1);
//
//         r_model_part1.CreateNewElement("Element3D8N", 1, {1, 2, 3, 4, 5, 6, 7, 8}, p_properties_0);
//         r_model_part2.CreateNewElement("Element3D8N", 2, {9, 10, 11, 12, 13, 14, 15, 16}, p_properties_1);
//
//         GidIO<> gid_io("test", GiD_PostBinary, SingleFile, WriteDeformed, WriteConditions);
//         gid_io.InitializeMesh(0.0);
//         gid_io.WriteMesh(r_main_model_part.GetMesh());
//         gid_io.FinalizeMesh();
//         gid_io.InitializeResults(0, r_main_model_part.GetMesh());
//         gid_io.FinalizeResults();

        has_intersection_reference = false;

        for (auto& r_point : first_hexa) {
            if (second_hexa.IsInside(r_point.Coordinates(), local_coords)) {
                has_intersection_reference = true;
                break;
            }
        }
        for (auto& r_point : second_hexa) {
            if (first_hexa.IsInside(r_point.Coordinates(), local_coords)) {
                has_intersection_reference = true;
                break;
            }
        }

        has_intersection = first_obb.HasIntersection(second_obb);

        KRATOS_EXPECT_EQ(has_intersection_reference, has_intersection);
    }

    /**
     *  Here we check the HasIntersection function of the OBB 3D
     * We compare with the intersection of two quadrilaterals
     */
    KRATOS_TEST_CASE_IN_SUITE(OBBHasIntersection3D, KratosCoreFastSuite)
    {
        array_1d<double, 3> first_center;
        first_center[0] = 0.5;
        first_center[1] = 0.5;
        first_center[2] = 0.1;

        array_1d<array_1d<double, 3>, 3> first_directions;
        first_directions[0][0] = std::sqrt(3.0)/3.0;
        first_directions[0][1] = std::sqrt(3.0)/3.0;
        first_directions[0][2] = std::sqrt(3.0)/3.0;
        MathUtils<double>::OrthonormalBasis(first_directions[0], first_directions[1], first_directions[2]);

        array_1d<double, 3> first_half_lenghts;
        first_half_lenghts[0] = std::sqrt(3.0)/3.0;
        first_half_lenghts[1] = 0.1;
        first_half_lenghts[2] = 0.1;

        OrientedBoundingBox<3> first_obb(first_center, first_directions, first_half_lenghts);

        array_1d<double, 3> second_center;
        second_center[0] = 0.0;
        second_center[1] = 0.0;
        second_center[2] = 0.0;

        array_1d<array_1d<double, 3>, 3> second_directions;
        second_directions[0][0] = 1.0;
        second_directions[0][1] = 0.0;
        second_directions[0][2] = 0.0;
        second_directions[1][0] = 0.0;
        second_directions[1][1] = 1.0;
        second_directions[1][2] = 0.0;
        second_directions[2][0] = 0.0;
        second_directions[2][1] = 0.0;
        second_directions[2][2] = 1.0;

        array_1d<double, 3> second_half_lenghts;
        second_half_lenghts[0] = 0.5;
        second_half_lenghts[1] = 0.1;
        second_half_lenghts[2] = 0.1;

        OrientedBoundingBox<3> second_obb(second_center, second_directions, second_half_lenghts);

        auto first_hexa = first_obb.GetEquivalentGeometry();
        auto second_hexa = second_obb.GetEquivalentGeometry();

        bool has_intersection_reference = false;

        array_1d<double, 3> local_coords;
        for (auto& r_point : first_hexa) {
            if (second_hexa.IsInside(r_point.Coordinates(), local_coords)) {
                has_intersection_reference = true;
                break;
            }
        }
        for (auto& r_point : second_hexa) {
            if (first_hexa.IsInside(r_point.Coordinates(), local_coords)) {
                has_intersection_reference = true;
                break;
            }
        }

//         // Debug
//         first_obb.GetEquivalentRotatedGeometry(second_hexa);
//         first_obb.GetEquivalentRotatedGeometry(first_hexa);
//
//         Model current_model;
//         ModelPart& r_main_model_part = current_model.CreateModelPart("Main");
//         ModelPart& r_model_part1 = r_main_model_part.CreateSubModelPart("1");
//         ModelPart& r_model_part2 = r_main_model_part.CreateSubModelPart("2");
//
//         for (int i = 0; i < 8; ++i) {
//             r_model_part1.CreateNewNode(i + 1, first_hexa[i].X(), first_hexa[i].Y(), first_hexa[i].Z());
//             r_model_part2.CreateNewNode(i + 9, second_hexa[i].X(), second_hexa[i].Y(), second_hexa[i].Z());
//         }
//
//         Properties::Pointer p_properties_0 = Kratos::make_shared<Properties>(0);
//         Properties::Pointer p_properties_1 = Kratos::make_shared<Properties>(1);
//
//         r_model_part1.CreateNewElement("Element3D8N", 1, {1, 2, 3, 4, 5, 6, 7, 8}, p_properties_0);
//         r_model_part2.CreateNewElement("Element3D8N", 2, {9, 10, 11, 12, 13, 14, 15, 16}, p_properties_1);
//
//         GidIO<> gid_io("test", GiD_PostBinary, SingleFile, WriteDeformed, WriteConditions);
//         gid_io.InitializeMesh(0.0);
//         gid_io.WriteMesh(r_main_model_part.GetMesh());
//         gid_io.FinalizeMesh();
//         gid_io.InitializeResults(0, r_main_model_part.GetMesh());
//         gid_io.FinalizeResults();

        bool has_intersection = first_obb.HasIntersection(second_obb);

        KRATOS_EXPECT_EQ(has_intersection_reference, has_intersection);

        // Moving up tyhe second OBB
        second_center[1] = 0.5;
        second_obb.SetCenter(second_center);
        second_hexa = second_obb.GetEquivalentGeometry();

//         // Debug
//         first_obb.GetEquivalentRotatedGeometry(second_hexa);
//         first_obb.GetEquivalentRotatedGeometry(first_hexa);
//
//         Model current_model;
//         ModelPart& r_main_model_part = current_model.CreateModelPart("Main");
//         ModelPart& r_model_part1 = r_main_model_part.CreateSubModelPart("1");
//         ModelPart& r_model_part2 = r_main_model_part.CreateSubModelPart("2");
//
//         for (int i = 0; i < 8; ++i) {
//             r_model_part1.CreateNewNode(i + 1, first_hexa[i].X(), first_hexa[i].Y(), first_hexa[i].Z());
//             r_model_part2.CreateNewNode(i + 9, second_hexa[i].X(), second_hexa[i].Y(), second_hexa[i].Z());
//         }
//
//         Properties::Pointer p_properties_0 = Kratos::make_shared<Properties>(0);
//         Properties::Pointer p_properties_1 = Kratos::make_shared<Properties>(1);
//
//         r_model_part1.CreateNewElement("Element3D8N", 1, {1, 2, 3, 4, 5, 6, 7, 8}, p_properties_0);
//         r_model_part2.CreateNewElement("Element3D8N", 2, {9, 10, 11, 12, 13, 14, 15, 16}, p_properties_1);
//
//         GidIO<> gid_io("test", GiD_PostBinary, SingleFile, WriteDeformed, WriteConditions);
//         gid_io.InitializeMesh(0.0);
//         gid_io.WriteMesh(r_main_model_part.GetMesh());
//         gid_io.FinalizeMesh();
//         gid_io.InitializeResults(0, r_main_model_part.GetMesh());
//         gid_io.FinalizeResults();

        has_intersection_reference = false;

        for (auto& r_point : first_hexa) {
            if (second_hexa.IsInside(r_point.Coordinates(), local_coords)) {
                has_intersection_reference = true;
                break;
            }
        }
        for (auto& r_point : second_hexa) {
            if (first_hexa.IsInside(r_point.Coordinates(), local_coords)) {
                has_intersection_reference = true;
                break;
            }
        }

        has_intersection = first_obb.HasIntersection(second_obb);

        KRATOS_EXPECT_EQ(has_intersection_reference, has_intersection);
    }

    /**
     *  Here we check the HasIntersection function of the OBB 3D
     * We compare with the intersection of two quadrilaterals
     */
    KRATOS_TEST_CASE_IN_SUITE(OBBHasIntersection3DComplexCase, KratosCoreFastSuite)
    {
        array_1d<double, 3> first_center;
        first_center[0] = 1.75;
        first_center[1] = 0;
        first_center[2] = 0.5;

        array_1d<array_1d<double, 3>, 3> first_directions;
        first_directions[0][0] = 0.447214;
        first_directions[0][1] = 0;
        first_directions[0][2] = 0.894427;
        first_directions[1][0] = 0.894427;
        first_directions[1][1] = -0;
        first_directions[1][2] = -0.447214;
        first_directions[2][0] = 0;
        first_directions[2][1] = 1;
        first_directions[2][2] = -0;

        array_1d<double, 3> first_half_lenghts;
        first_half_lenghts[0] = 0.659017;
        first_half_lenghts[1] = 0.547214;
        first_half_lenghts[2] = 0.1;

        OrientedBoundingBox<3> first_obb(first_center, first_directions, first_half_lenghts);

        array_1d<double, 3> second_center;
        second_center[0] = 1.74078;
        second_center[1] = 0.261432;
        second_center[2] = 0.5;

        array_1d<array_1d<double, 3>, 3> second_directions;
        second_directions[0][0] = 0.349731;
        second_directions[0][1] = 0.10609;
        second_directions[0][2] = 0.930824;
        second_directions[1][0] = 0.93685;
        second_directions[1][1] = -0.039604;
        second_directions[1][2] = -0.347482;
        second_directions[2][0] = 0;
        second_directions[2][1] = 0.993568;
        second_directions[2][2] = -0.113241;

        array_1d<double, 3> second_half_lenghts;
        second_half_lenghts[0] = 0.637159;
        second_half_lenghts[1] = 0.447482;
        second_half_lenghts[2] = 0.213241;

        OrientedBoundingBox<3> second_obb(second_center, second_directions, second_half_lenghts);

        auto first_hexa = first_obb.GetEquivalentGeometry();
        auto second_hexa = second_obb.GetEquivalentGeometry();

        bool has_intersection_reference = false;

        array_1d<double, 3> local_coords;
        for (auto& r_point : first_hexa) {
            if (second_hexa.IsInside(r_point.Coordinates(), local_coords)) {
                has_intersection_reference = true;
                break;
            }
        }
        for (auto& r_point : second_hexa) {
            if (first_hexa.IsInside(r_point.Coordinates(), local_coords)) {
                has_intersection_reference = true;
                break;
            }
        }

//         // Debug
//         first_obb.GetEquivalentRotatedGeometry(second_hexa);
//         first_obb.GetEquivalentRotatedGeometry(first_hexa);
//
//         Model current_model;
//         ModelPart& r_main_model_part = current_model.CreateModelPart("Main");
//         ModelPart& r_model_part1 = r_main_model_part.CreateSubModelPart("1");
//         ModelPart& r_model_part2 = r_main_model_part.CreateSubModelPart("2");
//
//         for (int i = 0; i < 8; ++i) {
//             r_model_part1.CreateNewNode(i + 1, first_hexa[i].X(), first_hexa[i].Y(), first_hexa[i].Z());
//             r_model_part2.CreateNewNode(i + 9, second_hexa[i].X(), second_hexa[i].Y(), second_hexa[i].Z());
//         }
//
//         Properties::Pointer p_properties_0 = Kratos::make_shared<Properties>(0);
//         Properties::Pointer p_properties_1 = Kratos::make_shared<Properties>(1);
//
//         r_model_part1.CreateNewElement("Element3D8N", 1, {1, 2, 3, 4, 5, 6, 7, 8}, p_properties_0);
//         r_model_part2.CreateNewElement("Element3D8N", 2, {9, 10, 11, 12, 13, 14, 15, 16}, p_properties_1);
//
//         GidIO<> gid_io("test", GiD_PostBinary, SingleFile, WriteDeformed, WriteConditions);
//         gid_io.InitializeMesh(0.0);
//         gid_io.WriteMesh(r_main_model_part.GetMesh());
//         gid_io.FinalizeMesh();
//         gid_io.InitializeResults(0, r_main_model_part.GetMesh());
//         gid_io.FinalizeResults();

        bool has_intersection = first_obb.HasIntersection(second_obb);

        KRATOS_EXPECT_EQ(has_intersection_reference, has_intersection);

        array_1d<double, 3> third_center;
        third_center[0] = 0.97685;
        third_center[1] = 0.0833096;
        third_center[2] = 0.5;

        array_1d<array_1d<double, 3>, 3> third_directions;
        third_directions[0][0] = 0.360587;
        third_directions[0][1] = 0.0595331;
        third_directions[0][2] = 0.930824;
        third_directions[1][0] = 0.932726;
        third_directions[1][1] = -0.0230152;
        third_directions[1][2] = -0.359852;
        third_directions[2][0] = 0;
        third_directions[2][1] = 0.997961;
        third_directions[2][2] = -0.0638271;

        array_1d<double, 3> third_half_lenghts;
        third_half_lenghts[0] = 0.637159;
        third_half_lenghts[1] = 0.459852;
        third_half_lenghts[2] = 0.163827;

        OrientedBoundingBox<3> third_obb(third_center, third_directions, third_half_lenghts);

        auto third_hexa = third_obb.GetEquivalentGeometry();

        has_intersection_reference = false;

        for (auto& r_point : first_hexa) {
            if (third_hexa.IsInside(r_point.Coordinates(), local_coords)) {
                has_intersection_reference = true;
                break;
            }
        }
        for (auto& r_point : third_hexa) {
            if (first_hexa.IsInside(r_point.Coordinates(), local_coords)) {
                has_intersection_reference = true;
                break;
            }
        }

        has_intersection = first_obb.HasIntersection(third_obb);

        KRATOS_EXPECT_EQ(has_intersection_reference, has_intersection);
    }

    }
}  // namespace Kratos.
