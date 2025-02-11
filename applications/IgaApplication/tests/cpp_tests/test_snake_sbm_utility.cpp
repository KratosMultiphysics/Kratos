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
#include "utilities/nurbs_utilities/snake_sbm_utilities.h"
#include "includes/kratos_parameters.h"

namespace Kratos::Testing
{

// Tests the SnakeSbmUtilities with a square outer geometry
KRATOS_TEST_CASE_IN_SUITE(SnakeSbmUtilitySquareOuter, KratosIgaFastSuite)
{
    
    Model model;
    ModelPart& iga_model_part = model.CreateModelPart("iga_model_part");

    iga_model_part.CreateSubModelPart("surrogate_inner");
    ModelPart& surrogate_sub_model_part_outer = iga_model_part.CreateSubModelPart("surrogate_outer");

    ModelPart& skin_model_part_inner_initial = model.CreateModelPart("skin_model_part_inner_initial");
    ModelPart& skin_model_part_outer_initial = model.CreateModelPart("skin_model_part_outer_initial");

    skin_model_part_outer_initial.CreateNewProperties(0);

    skin_model_part_outer_initial.CreateNewNode(1, 0.0, 0.0, 0.0);
    skin_model_part_outer_initial.CreateNewNode(2, 2.0, 0.0, 0.0);
    skin_model_part_outer_initial.CreateNewNode(3, 2.0, 2.0, 0.0);
    skin_model_part_outer_initial.CreateNewNode(4, 0.0, 2.0, 0.0);
        
    Properties::Pointer p_prop = skin_model_part_outer_initial.pGetProperties(0);
    skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 1, {{1, 2}}, p_prop);
    skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 2, {{2, 3}}, p_prop);
    skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 3, {{3, 4}}, p_prop);
    skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 4, {{4, 1}}, p_prop);

    ModelPart& skin_model_part = model.CreateModelPart("skin_model_part");
    skin_model_part.CreateSubModelPart("inner");
    skin_model_part.CreateSubModelPart("outer");
    
    const std::vector<double> list_knot_u = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0};
    const std::vector<double> list_knot_v = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0};

    Vector unique_knot_vector_u(list_knot_u.size());
    Vector unique_knot_vector_v(list_knot_v.size());

    // Fill the vectors using a loop
    for (std::size_t i = 0; i < list_knot_u.size(); ++i) 
        unique_knot_vector_u[i] = list_knot_u[i];
    for (std::size_t i = 0; i < list_knot_v.size(); ++i) 
        unique_knot_vector_v[i] = list_knot_v[i];

    Kratos::Parameters mParameters(R"(
        {
            "sbm_parameters": {
                "lambda_outer": 0.5,
                "number_of_inner_loops": 0
            }
        }
    )");
        
    SnakeSbmUtilities::CreateTheSnakeCoordinates(iga_model_part, skin_model_part_inner_initial, skin_model_part_outer_initial, skin_model_part, 0,
        unique_knot_vector_u, unique_knot_vector_v, mParameters) ;
    
    const double tolerance = 1.0e-6;

    // Expected coordinates of nodes (modify according to actual expected values)
    std::vector<std::array<double, 3>> expected_coordinates = {
        {0, 0, 0},
        {2, 0, 0},
        {0, 0.4,0},
        {2, 0.4,0},
        {0, 0.8,0},
        {2, 0.8,0},
        {0, 1.2,0},
        {2, 1.2,0},
        {0, 1.6,0},
        {2, 1.6,0},
        {0, 0, 0},
        {0, 2, 0},
        {0.4, 0, 0},
        {0.4, 2, 0},
        {0.8, 0, 0},
        {0.8, 2, 0},
        {1.2, 0, 0},
        {1.2, 2, 0},
        {1.6, 0, 0},
        {1.6, 2, 0}
    };
    
    // Ensure the number of nodes matches expectation
    KRATOS_EXPECT_NEAR(surrogate_sub_model_part_outer.NumberOfConditions(), expected_coordinates.size(), tolerance);
    
    // Iterate over nodes and compare coordinates
    unsigned int i = 0;
    for (auto& cond : surrogate_sub_model_part_outer.Conditions()) {
        std::array<double, 3> node_coords = {cond.GetGeometry()[0].X(), cond.GetGeometry()[0].Y(), cond.GetGeometry()[0].Z()};
        // Compare each coordinate
        KRATOS_EXPECT_NEAR(expected_coordinates[i][0], node_coords[0], tolerance);
        KRATOS_EXPECT_NEAR(expected_coordinates[i][1], node_coords[1], tolerance);
        KRATOS_EXPECT_NEAR(expected_coordinates[i][2], node_coords[2], tolerance);
        i++;
    }

}

// Tests the SnakeSbmUtilities with an inner geometry
KRATOS_TEST_CASE_IN_SUITE(SnakeSbmUtilityInner, KratosIgaFastSuite)
{
    
    Model model;
    ModelPart& iga_model_part = model.CreateModelPart("iga_model_part");

    ModelPart& surrogate_sub_model_part_inner = iga_model_part.CreateSubModelPart("surrogate_inner");
    iga_model_part.CreateSubModelPart("surrogate_outer");

    ModelPart& skin_model_part_inner_initial = model.CreateModelPart("skin_model_part_inner_initial");
    ModelPart& skin_model_part_outer_initial = model.CreateModelPart("skin_model_part_outer_initial");

    skin_model_part_inner_initial.CreateNewProperties(0);

    // object 1
    skin_model_part_inner_initial.CreateNewNode(1, 0.5, 0.5, 0.0);
    skin_model_part_inner_initial.CreateNewNode(2, 3.5, 0.5, 0.0);
    skin_model_part_inner_initial.CreateNewNode(3, 3.5, 1.0, 0.0);
    skin_model_part_inner_initial.CreateNewNode(4, 3.0, 1.0, 0.0);
    skin_model_part_inner_initial.CreateNewNode(5, 3.0, 2.0, 0.0);
    skin_model_part_inner_initial.CreateNewNode(6, 3.5, 2.0, 0.0);
    skin_model_part_inner_initial.CreateNewNode(7, 3.5, 3.5, 0.0);
    skin_model_part_inner_initial.CreateNewNode(8, 3.0, 3.5, 0.0);
    skin_model_part_inner_initial.CreateNewNode(9, 3.0, 3.0, 0.0);
    skin_model_part_inner_initial.CreateNewNode(10, 2.0, 3.0, 0.0);
    skin_model_part_inner_initial.CreateNewNode(11, 2.0, 3.5, 0.0);
    skin_model_part_inner_initial.CreateNewNode(12, 0.5, 3.5, 0.0);

    Properties::Pointer p_prop = skin_model_part_inner_initial.pGetProperties(0);
    skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 1, {{1, 2}}, p_prop);
    skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 2, {{2, 3}}, p_prop);
    skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 3, {{3, 4}}, p_prop);
    skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 4, {{4, 5}}, p_prop);
    skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 5, {{5, 6}}, p_prop);
    skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 6, {{6, 7}}, p_prop);
    skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 7, {{7, 8}}, p_prop);
    skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 8, {{8, 9}}, p_prop);
    skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 9, {{9, 10}}, p_prop);
    skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 10, {{10, 11}}, p_prop);
    skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 11, {{11, 12}}, p_prop);
    skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 12, {{12, 1}}, p_prop);


    // object 2
    skin_model_part_inner_initial.CreateNewNode(13, 1.5, 4.6, 0.0);
    skin_model_part_inner_initial.CreateNewNode(14, 3.7, 4.7, 0.0);
    skin_model_part_inner_initial.CreateNewNode(15, 3.0, 5.5, 0.0);
    skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 13, {{13, 14}}, p_prop);
    skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 14, {{14, 15}}, p_prop);
    skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 15, {{15, 13}}, p_prop);

    ModelPart& skin_model_part = model.CreateModelPart("skin_model_part");
    skin_model_part.CreateSubModelPart("inner");
    skin_model_part.CreateSubModelPart("outer");
    
    const std::vector<double> list_knot_u = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0};
    const std::vector<double> list_knot_v = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.4, 4.8, 5.2, 5.6, 6.0};

    Vector unique_knot_vector_u(list_knot_u.size());
    Vector unique_knot_vector_v(list_knot_v.size());

    // Fill the vectors using a loop
    for (std::size_t i = 0; i < list_knot_u.size(); ++i) 
        unique_knot_vector_u[i] = list_knot_u[i];
    for (std::size_t i = 0; i < list_knot_v.size(); ++i) 
        unique_knot_vector_v[i] = list_knot_v[i];

    Kratos::Parameters mParameters(R"(
        {
            "sbm_parameters": {
                    "lambda_inner": 0.5,
                    "number_of_inner_loops": 2
            }
        }
    )");

    SnakeSbmUtilities::CreateTheSnakeCoordinates(iga_model_part, skin_model_part_inner_initial, skin_model_part_outer_initial, skin_model_part, 0,
        unique_knot_vector_u, unique_knot_vector_v, mParameters) ;
    
    const double tolerance = 1.0e-6;

    // Expected coordinates of nodes (modify according to actual expected values)
    std::vector<std::array<double, 3>> expected_coordinates = {
            {0.4, 0.4, 0},
            {3.6, 0.4, 0},
            {0.4, 0.8, 0},
            {3.2, 0.8, 0},
            {0.4, 1.2, 0},
            {3.2, 1.2, 0},
            {0.4, 1.6, 0},
            {3.2, 1.6, 0},
            {0.4, 2, 0},
            {3.6, 2, 0},
            {0.4, 2.4, 0},
            {3.6, 2.4, 0},
            {0.4, 2.8, 0},
            {3.6, 2.8, 0},
            {0.4, 3.2, 0},
            {2, 3.2, 0},
            {3.2, 3.2, 0},
            {3.6, 3.2, 0},
            {0.4, 0.4, 0},
            {0.4, 3.6, 0},
            {0.8, 0.4, 0},
            {0.8, 3.6, 0},
            {1.2, 0.4, 0},
            {1.2, 3.6, 0},
            {1.6, 0.4, 0},
            {1.6, 3.6, 0},
            {2, 0.4, 0},
            {2, 3.2, 0},
            {2.4, 0.4, 0},
            {2.4, 3.2, 0},
            {2.8, 0.4, 0},
            {2.8, 3.2, 0},
            {3.2, 0.4, 0},
            {3.2, 0.8, 0},
            {3.2, 2, 0},
            {3.2, 3.6, 0},
            {2, 4.8, 0},
            {3.6, 4.8, 0},
            {2.8, 5.2, 0},
            {3.2, 5.2, 0},
            {2, 4.8, 0},
            {2, 5.2, 0},
            {2.4, 4.8, 0},
            {2.4, 5.2, 0},
            {2.8, 4.8, 0},
            {2.8, 5.6, 0},
            {3.2, 4.8, 0},
            {3.2, 5.2, 0}
    };
    
    // Ensure the number of nodes matches expectation
    KRATOS_EXPECT_NEAR(surrogate_sub_model_part_inner.NumberOfConditions(), expected_coordinates.size(), tolerance);
    
    // Iterate over nodes and compare coordinates
    unsigned int i = 0;
    for (auto& cond : surrogate_sub_model_part_inner.Conditions()) {
        std::array<double, 3> node_coords = {cond.GetGeometry()[0].X(), cond.GetGeometry()[0].Y(), cond.GetGeometry()[0].Z()};
        // Compare each coordinate
        KRATOS_EXPECT_NEAR(expected_coordinates[i][0], node_coords[0], tolerance);
        KRATOS_EXPECT_NEAR(expected_coordinates[i][1], node_coords[1], tolerance);
        KRATOS_EXPECT_NEAR(expected_coordinates[i][2], node_coords[2], tolerance);
        i++;
    }
}

// Tests the SnakeSbmUtilities with an inner and an outer geometry
KRATOS_TEST_CASE_IN_SUITE(SnakeSbmUtilityInnerOuter, KratosIgaFastSuite)
{
    
    Model model;
    ModelPart& iga_model_part = model.CreateModelPart("iga_model_part");

    ModelPart& surrogate_sub_model_part_inner = iga_model_part.CreateSubModelPart("surrogate_inner");
    ModelPart& surrogate_sub_model_part_outer = iga_model_part.CreateSubModelPart("surrogate_outer");

    ModelPart& skin_model_part_inner_initial = model.CreateModelPart("skin_model_part_inner_initial");
    ModelPart& skin_model_part_outer_initial = model.CreateModelPart("skin_model_part_outer_initial");

    skin_model_part_inner_initial.CreateNewProperties(0);
    skin_model_part_outer_initial.CreateNewProperties(0);
    
    // inner
    skin_model_part_inner_initial.CreateNewNode(1, 1.5, 0.5, 0.0);
    skin_model_part_inner_initial.CreateNewNode(2, 3.0, 0.5, 0.0);
    skin_model_part_inner_initial.CreateNewNode(3, 3.0, 3.5, 0.0);
    skin_model_part_inner_initial.CreateNewNode(4, 1.5, 3.5, 0.0);
    Properties::Pointer p_prop = skin_model_part_inner_initial.pGetProperties(0);
    skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 1, {{1, 2}}, p_prop);
    skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 2, {{2, 3}}, p_prop);
    skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 3, {{3, 4}}, p_prop);
    skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 4, {{4, 1}}, p_prop);

    // outer
    p_prop = skin_model_part_outer_initial.pGetProperties(0);
    skin_model_part_outer_initial.CreateNewNode(1, 0.4, 0.0, 0.0);
    skin_model_part_outer_initial.CreateNewNode(2, 3.5, 0.0, 0.0);
    skin_model_part_outer_initial.CreateNewNode(3, 3.2, 5.0, 0.0);
    skin_model_part_outer_initial.CreateNewNode(4, 2.5, 5.0, 0.0);
    skin_model_part_outer_initial.CreateNewNode(5, 2.5, 4.5, 0.0);
    skin_model_part_outer_initial.CreateNewNode(6, 1.5, 4.5, 0.0);
    skin_model_part_outer_initial.CreateNewNode(7, 1.1, 4.9, 0.0);
    skin_model_part_outer_initial.CreateNewNode(8, 1.0, 5.5, 0.0);
    skin_model_part_outer_initial.CreateNewNode(9, 0.9, 5.1, 0.0);
    skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 1, {{1, 2}}, p_prop);
    skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 2, {{2, 3}}, p_prop);
    skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 3, {{3, 4}}, p_prop);
    skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 4, {{4, 5}}, p_prop);
    skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 5, {{5, 6}}, p_prop);
    skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 6, {{6, 7}}, p_prop);
    skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 7, {{7, 8}}, p_prop);
    skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 8, {{8, 9}}, p_prop);
    skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 9, {{9, 1}}, p_prop);


    ModelPart& skin_model_part = model.CreateModelPart("skin_model_part");
    skin_model_part.CreateSubModelPart("inner");
    skin_model_part.CreateSubModelPart("outer");
    
    const std::vector<double> list_knot_u = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0};
    const std::vector<double> list_knot_v = {0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.3, 3.6, 3.9, 4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0};

    Vector unique_knot_vector_u(list_knot_u.size());
    Vector unique_knot_vector_v(list_knot_v.size());

    // Fill the vectors using a loop
    for (std::size_t i = 0; i < list_knot_u.size(); ++i) 
        unique_knot_vector_u[i] = list_knot_u[i];
    for (std::size_t i = 0; i < list_knot_v.size(); ++i) 
        unique_knot_vector_v[i] = list_knot_v[i];

    Kratos::Parameters mParameters(R"(
        {
            "sbm_parameters": {
                    "lambda_inner": 1.0,
                    "lambda_outer": 0.5,
                    "number_of_inner_loops": 1
            }
        }
    )");
        
    SnakeSbmUtilities::CreateTheSnakeCoordinates(iga_model_part, skin_model_part_inner_initial, skin_model_part_outer_initial, skin_model_part, 0,
        unique_knot_vector_u, unique_knot_vector_v, mParameters) ;
    
    const double tolerance = 1.0e-6;

    std::vector<std::array<double, 3>> expected_coordinates_outer = { 
        {0.4, 0, 0},
        {3.6, 0, 0},
        {0.4, 0.3, 0},
        {3.6, 0.3, 0},
        {0.4, 0.6, 0},
        {3.6, 0.6, 0},
        {0.4, 0.9, 0},
        {3.6, 0.9, 0},
        {0.4, 1.2, 0},
        {3.6, 1.2, 0},
        {0.4, 1.5, 0},
        {3.6, 1.5, 0},
        {0.4, 1.8, 0},
        {3.2, 1.8, 0},
        {0.8, 2.1, 0},
        {3.2, 2.1, 0},
        {0.8, 2.4, 0},
        {3.2, 2.4, 0},
        {0.8, 2.7, 0},
        {3.2, 2.7, 0},
        {0.8, 3, 0},
        {3.2, 3, 0},
        {0.8, 3.3, 0},
        {3.2, 3.3, 0},
        {0.8, 3.6, 0},
        {3.2, 3.6, 0},
        {0.8, 3.9, 0},
        {3.2, 3.9, 0},
        {0.8, 4.2, 0},
        {3.2, 4.2, 0},
        {0.8, 4.5, 0},
        {1.2, 4.5, 0},
        {2.4, 4.5, 0},
        {3.2, 4.5, 0},
        {0.8, 4.8, 0},
        {1.2, 4.8, 0},
        {2.4, 4.8, 0},
        {3.2, 4.8, 0},
        {0.4, 0, 0},
        {0.4, 2.1, 0},
        {0.8, 0, 0},
        {0.8, 5.1, 0},
        {1.2, 0, 0},
        {1.2, 4.5, 0},
        {1.6, 0, 0},
        {1.6, 4.5, 0},
        {2, 0, 0},
        {2, 4.5, 0},
        {2.4, 0, 0},
        {2.4, 5.1, 0},
        {2.8, 0, 0},
        {2.8, 5.1, 0},
        {3.2, 0, 0},
        {3.2, 1.8, 0},
    };
    
    std::vector<std::array<double, 3>> expected_coordinates_inner = {
        {1.6, 0.6, 0},
        {2.8, 0.6, 0},
        {1.6, 0.9, 0},
        {2.8, 0.9, 0},
        {1.6, 1.2, 0},
        {2.8, 1.2, 0},
        {1.6, 1.5, 0},
        {2.8, 1.5, 0},
        {1.6, 1.8, 0},
        {2.8, 1.8, 0},
        {1.6, 2.1, 0},
        {2.8, 2.1, 0},
        {1.6, 2.4, 0},
        {2.8, 2.4, 0},
        {1.6, 2.7, 0},
        {2.8, 2.7, 0},
        {1.6, 3, 0},
        {2.8, 3, 0},
        {1.6, 0.6, 0},
        {1.6, 3.3, 0},
        {2, 0.6, 0},
        {2, 3.3, 0},
        {2.4, 0.6, 0},
        {2.4, 3.3, 0}
    };
    
    // Ensure the number of nodes matches expectation
    KRATOS_EXPECT_NEAR(surrogate_sub_model_part_inner.NumberOfConditions(), expected_coordinates_inner.size(), tolerance);
    KRATOS_EXPECT_NEAR(surrogate_sub_model_part_outer.NumberOfConditions(), expected_coordinates_outer.size(), tolerance);
    
    // Iterate over nodes and compare coordinates
    unsigned int i = 0;
    for (auto& cond : surrogate_sub_model_part_inner.Conditions()) {
        std::array<double, 3> node_coords = {cond.GetGeometry()[0].X(), cond.GetGeometry()[0].Y(), cond.GetGeometry()[0].Z()};
        // Compare each coordinate
        KRATOS_EXPECT_NEAR(expected_coordinates_inner[i][0], node_coords[0], tolerance);
        KRATOS_EXPECT_NEAR(expected_coordinates_inner[i][1], node_coords[1], tolerance);
        KRATOS_EXPECT_NEAR(expected_coordinates_inner[i][2], node_coords[2], tolerance);
        i++;
    }
    i = 0;
    for (auto& cond : surrogate_sub_model_part_outer.Conditions()) {
        std::array<double, 3> node_coords = {cond.GetGeometry()[0].X(), cond.GetGeometry()[0].Y(), cond.GetGeometry()[0].Z()};
        // Compare each coordinate
        KRATOS_EXPECT_NEAR(expected_coordinates_outer[i][0], node_coords[0], tolerance);
        KRATOS_EXPECT_NEAR(expected_coordinates_outer[i][1], node_coords[1], tolerance);
        KRATOS_EXPECT_NEAR(expected_coordinates_outer[i][2], node_coords[2], tolerance);
        i++;
    }

}

}
