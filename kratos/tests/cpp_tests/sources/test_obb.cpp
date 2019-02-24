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
#include "includes/obb.h"
#include "includes/checks.h"
#include "geometries/point.h"

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

        OBB<2> first_obb(first_center, first_directions, first_half_lenghts);

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

        OBB<2> second_obb(second_center, second_directions, second_half_lenghts);

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

        KRATOS_CHECK_EQUAL(has_intersection_reference, has_intersection);

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

        KRATOS_CHECK_EQUAL(has_intersection_reference, has_intersection);
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

        OBB<3> first_obb(first_center, first_directions, first_half_lenghts);

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

        OBB<3> second_obb(second_center, second_directions, second_half_lenghts);

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

        KRATOS_CHECK_EQUAL(has_intersection_reference, has_intersection);

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

        KRATOS_CHECK_EQUAL(has_intersection_reference, has_intersection);
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

        OBB<3> first_obb(first_center, first_directions, first_half_lenghts);

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

        OBB<3> second_obb(second_center, second_directions, second_half_lenghts);

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

        KRATOS_CHECK_EQUAL(has_intersection_reference, has_intersection);

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

        KRATOS_CHECK_EQUAL(has_intersection_reference, has_intersection);
    }

    }
}  // namespace Kratos.
