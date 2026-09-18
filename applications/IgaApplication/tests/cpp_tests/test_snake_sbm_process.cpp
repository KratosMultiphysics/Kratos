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
#include "custom_processes/snake_sbm_process.h"
#include "includes/kratos_parameters.h"
#include "iga_application_variables.h"
#include "geometries/nurbs_curve_geometry.h"

namespace Kratos::Testing
{

// Tests the SnakeSbmUProcess with a square outer geometry
KRATOS_TEST_CASE_IN_SUITE(SnakeSbmProcessSquareOuter, KratosIgaFastSuite)
{
    
    Model model;
    ModelPart& iga_model_part = model.CreateModelPart("iga_model_part");

    iga_model_part.CreateSubModelPart("surrogate_inner");
    ModelPart& surrogate_sub_model_part_outer = iga_model_part.CreateSubModelPart("surrogate_outer");

    model.CreateModelPart("skin_model_part_inner_initial");
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

    Kratos::Parameters snake_parameters(R"(
        {
            "model_part_name" : "iga_model_part",
            "skin_model_part_inner_initial_name" : "skin_model_part_inner_initial",
            "skin_model_part_outer_initial_name" : "skin_model_part_outer_initial",
            "skin_model_part_name" : "skin_model_part",
            "echo_level" : 4,
            "lambda_inner" : 0.5,
            "lambda_outer" : 0.5,
            "number_of_inner_loops": 0
        }
    )");

    SnakeSbmProcess snake_sbm_process(model, snake_parameters);

    iga_model_part.SetValue(KNOT_VECTOR_U, unique_knot_vector_u);
    iga_model_part.SetValue(KNOT_VECTOR_V, unique_knot_vector_v);
    
    snake_sbm_process.Execute() ;
    
    const double tolerance = 1.0e-6;

    // Expected coordinates of nodes (modify according to actual expected values)
    std::vector<std::array<double, 3>> expected_coordinates = {
        {0, 0, 0},
        {2, 0, 0},
        {0, 0.4, 0},
        {2, 0.4, 0},
        {0, 0.8, 0},
        {2, 0.8, 0},
        {0, 1.2, 0},
        {2, 1.2, 0},
        {0, 1.6, 0},
        {2, 1.6, 0},
        {0.4, 0, 0},
        {0.4, 2, 0},
        {0.8, 0, 0},
        {0.8, 2, 0},
        {1.2, 0, 0},
        {1.2, 2, 0},
        {1.6, 0, 0},
        {1.6, 2, 0},
        {2, 0, 0},
        {2, 2, 0}
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

// Tests the SnakeSbmProcess with an inner geometry
KRATOS_TEST_CASE_IN_SUITE(SnakeSbmProcessInner, KratosIgaFastSuite)
{
    
    Model model;
    ModelPart& iga_model_part = model.CreateModelPart("iga_model_part");

    ModelPart& surrogate_sub_model_part_inner = iga_model_part.CreateSubModelPart("surrogate_inner");
    iga_model_part.CreateSubModelPart("surrogate_outer");

    ModelPart& skin_model_part_inner_initial = model.CreateModelPart("skin_model_part_inner_initial");
    model.CreateModelPart("skin_model_part_outer_initial");

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

    Kratos::Parameters snake_parameters(R"(
        {
            "model_part_name" : "iga_model_part",
            "skin_model_part_inner_initial_name" : "skin_model_part_inner_initial",
            "skin_model_part_outer_initial_name" : "skin_model_part_outer_initial",
            "skin_model_part_name" : "skin_model_part",
            "echo_level" : 0,
            "lambda_inner" : 0.5,
            "lambda_outer" : 0.5,
            "number_of_inner_loops": 2
        }
    )");

    SnakeSbmProcess snake_sbm_process(model, snake_parameters);

    iga_model_part.SetValue(KNOT_VECTOR_U, unique_knot_vector_u);
    iga_model_part.SetValue(KNOT_VECTOR_V, unique_knot_vector_v);

    snake_sbm_process.Execute() ;
    
    const double tolerance = 1.0e-6;

    // Expected coordinates of nodes (modify according to actual expected values)
    std::vector<std::array<double, 3>> expected_coordinates = {
            {0.4, 0.8, 0},
            {3.6, 0.8, 0},
            {0.4, 1.2, 0},
            {3.2, 1.2, 0},
            {0.4, 1.6, 0},
            {3.2, 1.6, 0},
            {0.4, 2.0, 0},
            {3.2, 2.0, 0},
            {0.4, 2.4, 0},
            {3.6, 2.4, 0},
            {0.4, 2.8, 0},
            {3.6, 2.8, 0},
            {0.4, 3.2, 0},
            {3.6, 3.2, 0},
            {0.4, 3.6, 0},
            {2, 3.6, 0},
            {3.2, 3.6, 0},
            {3.6, 3.6, 0},
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
            {2, 5.2, 0},
            {3.6, 5.2, 0},
            {2.8, 5.6, 0},
            {3.2, 5.6, 0},
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

// // Tests the SnakeSbmProcess with an inner and an outer geometry
// KRATOS_TEST_CASE_IN_SUITE(SnakeSbmProcessInnerOuter, KratosIgaFastSuite)
// {
    
//     Model model;
//     ModelPart& iga_model_part = model.CreateModelPart("iga_model_part");

//     ModelPart& surrogate_sub_model_part_inner = iga_model_part.CreateSubModelPart("surrogate_inner");
//     ModelPart& surrogate_sub_model_part_outer = iga_model_part.CreateSubModelPart("surrogate_outer");

//     ModelPart& skin_model_part_inner_initial = model.CreateModelPart("skin_model_part_inner_initial");
//     ModelPart& skin_model_part_outer_initial = model.CreateModelPart("skin_model_part_outer_initial");

//     skin_model_part_inner_initial.CreateNewProperties(0);
//     skin_model_part_outer_initial.CreateNewProperties(0);
    
//     // inner
//     skin_model_part_inner_initial.CreateNewNode(1, 1.5, 0.5, 0.0);
//     skin_model_part_inner_initial.CreateNewNode(2, 3.0, 0.5, 0.0);
//     skin_model_part_inner_initial.CreateNewNode(3, 3.0, 3.5, 0.0);
//     skin_model_part_inner_initial.CreateNewNode(4, 1.5, 3.5, 0.0);
//     Properties::Pointer p_prop = skin_model_part_inner_initial.pGetProperties(0);
//     skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 1, {{1, 2}}, p_prop);
//     skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 2, {{2, 3}}, p_prop);
//     skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 3, {{3, 4}}, p_prop);
//     skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 4, {{4, 1}}, p_prop);

//     // outer
//     p_prop = skin_model_part_outer_initial.pGetProperties(0);
//     skin_model_part_outer_initial.CreateNewNode(1, 0.4, 0.0, 0.0);
//     skin_model_part_outer_initial.CreateNewNode(2, 3.5, 0.0, 0.0);
//     skin_model_part_outer_initial.CreateNewNode(3, 3.2, 5.0, 0.0);
//     skin_model_part_outer_initial.CreateNewNode(4, 2.5, 5.0, 0.0);
//     skin_model_part_outer_initial.CreateNewNode(5, 2.5, 4.5, 0.0);
//     skin_model_part_outer_initial.CreateNewNode(6, 1.5, 4.5, 0.0);
//     skin_model_part_outer_initial.CreateNewNode(7, 1.1, 4.9, 0.0);
//     skin_model_part_outer_initial.CreateNewNode(8, 1.0, 5.5, 0.0);
//     skin_model_part_outer_initial.CreateNewNode(9, 0.9, 5.1, 0.0);
//     skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 1, {{1, 2}}, p_prop);
//     skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 2, {{2, 3}}, p_prop);
//     skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 3, {{3, 4}}, p_prop);
//     skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 4, {{4, 5}}, p_prop);
//     skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 5, {{5, 6}}, p_prop);
//     skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 6, {{6, 7}}, p_prop);
//     skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 7, {{7, 8}}, p_prop);
//     skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 8, {{8, 9}}, p_prop);
//     skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 9, {{9, 1}}, p_prop);


//     ModelPart& skin_model_part = model.CreateModelPart("skin_model_part");
//     skin_model_part.CreateSubModelPart("inner");
//     skin_model_part.CreateSubModelPart("outer");
    
//     const std::vector<double> list_knot_u = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0};
//     const std::vector<double> list_knot_v = {0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.3, 3.6, 3.9, 4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0};

//     Vector unique_knot_vector_u(list_knot_u.size());
//     Vector unique_knot_vector_v(list_knot_v.size());

//     // Fill the vectors using a loop
//     for (std::size_t i = 0; i < list_knot_u.size(); ++i) 
//         unique_knot_vector_u[i] = list_knot_u[i];
//     for (std::size_t i = 0; i < list_knot_v.size(); ++i) 
//         unique_knot_vector_v[i] = list_knot_v[i];

//     Kratos::Parameters snake_parameters(R"(
//         {
//             "model_part_name" : "iga_model_part",
//             "skin_model_part_inner_initial_name" : "skin_model_part_inner_initial",
//             "skin_model_part_outer_initial_name" : "skin_model_part_outer_initial",
//             "skin_model_part_name" : "skin_model_part",
//             "echo_level" : 0,
//             "lambda_inner" : 1.0,
//             "lambda_outer" : 0.5,
//             "number_of_inner_loops": 1
//         }
//     )");

//     SnakeSbmProcess snake_sbm_process(model, snake_parameters);

//     iga_model_part.SetValue(KNOT_VECTOR_U, unique_knot_vector_u);
//     iga_model_part.SetValue(KNOT_VECTOR_V, unique_knot_vector_v);

//     snake_sbm_process.Execute() ;

//     const double tolerance = 1.0e-6;

//     std::vector<std::array<double, 3>> expected_coordinates_outer = { 
//         {0.4, 0, 0},
//         {3.6, 0, 0},
//         {0.4, 0.3, 0},
//         {3.6, 0.3, 0},
//         {0.4, 0.6, 0},
//         {3.6, 0.6, 0},
//         {0.4, 0.9, 0},
//         {3.6, 0.9, 0},
//         {0.4, 1.2, 0},
//         {3.6, 1.2, 0},
//         {0.4, 1.5, 0},
//         {3.6, 1.5, 0},
//         {0.4, 1.8, 0},
//         {3.2, 1.8, 0},
//         {0.8, 2.1, 0},
//         {3.2, 2.1, 0},
//         {0.8, 2.4, 0},
//         {3.2, 2.4, 0},
//         {0.8, 2.7, 0},
//         {3.2, 2.7, 0},
//         {0.8, 3, 0},
//         {3.2, 3, 0},
//         {0.8, 3.3, 0},
//         {3.2, 3.3, 0},
//         {0.8, 3.6, 0},
//         {3.2, 3.6, 0},
//         {0.8, 3.9, 0},
//         {3.2, 3.9, 0},
//         {0.8, 4.2, 0},
//         {3.2, 4.2, 0},
//         {0.8, 4.5, 0},
//         {1.2, 4.5, 0},
//         {2.4, 4.5, 0},
//         {3.2, 4.5, 0},
//         {2.4, 4.8, 0},
//         {3.2, 4.8, 0},
//         {0.8, 0, 0},
//         {0.8, 2.1, 0},
//         {1.2, 0, 0},
//         {1.2, 4.8, 0},
//         {1.6, 0, 0},
//         {1.6, 4.5, 0},
//         {2, 0, 0},
//         {2, 4.5, 0},
//         {2.4, 0, 0},
//         {2.4, 4.5, 0},
//         {2.8, 0, 0},
//         {2.8, 5.1, 0},
//         {3.2, 0, 0},
//         {3.2, 5.1, 0},
//         {3.6, 0, 0},
//         {3.6, 1.8, 0}
//     };
    
//     std::vector<std::array<double, 3>> expected_coordinates_inner = {
//         {1.6, 0.9, 0},
//         {2.8, 0.9, 0},
//         {1.6, 1.2, 0},
//         {2.8, 1.2, 0},
//         {1.6, 1.5, 0},
//         {2.8, 1.5, 0},
//         {1.6, 1.8, 0},
//         {2.8, 1.8, 0},
//         {1.6, 2.1, 0},
//         {2.8, 2.1, 0},
//         {1.6, 2.4, 0},
//         {2.8, 2.4, 0},
//         {1.6, 2.7, 0},
//         {2.8, 2.7, 0},
//         {1.6, 3, 0},
//         {2.8, 3, 0},
//         {1.6, 3.3, 0},
//         {2.8, 3.3, 0},
//         {1.6, 0.6, 0},
//         {1.6, 3.3, 0},
//         {2, 0.6, 0},
//         {2, 3.3, 0},
//         {2.4, 0.6, 0},
//         {2.4, 3.3, 0}
//     };
    
//     // Ensure the number of nodes matches expectation
//     KRATOS_EXPECT_NEAR(surrogate_sub_model_part_inner.NumberOfConditions(), expected_coordinates_inner.size(), tolerance);
//     KRATOS_EXPECT_NEAR(surrogate_sub_model_part_outer.NumberOfConditions(), expected_coordinates_outer.size(), tolerance);
    
//     // Iterate over nodes and compare coordinates
//     unsigned int i = 0;
//     for (auto& cond : surrogate_sub_model_part_inner.Conditions()) {
//         std::array<double, 3> node_coords = {cond.GetGeometry()[0].X(), cond.GetGeometry()[0].Y(), cond.GetGeometry()[0].Z()};
//         // Compare each coordinate
//         KRATOS_EXPECT_NEAR(expected_coordinates_inner[i][0], node_coords[0], tolerance);
//         KRATOS_EXPECT_NEAR(expected_coordinates_inner[i][1], node_coords[1], tolerance);
//         KRATOS_EXPECT_NEAR(expected_coordinates_inner[i][2], node_coords[2], tolerance);
//         i++;
//     }
//     i = 0;
//     for (auto& cond : surrogate_sub_model_part_outer.Conditions()) {
//         std::array<double, 3> node_coords = {cond.GetGeometry()[0].X(), cond.GetGeometry()[0].Y(), cond.GetGeometry()[0].Z()};
//         // Compare each coordinate
//         KRATOS_EXPECT_NEAR(expected_coordinates_outer[i][0], node_coords[0], tolerance);
//         KRATOS_EXPECT_NEAR(expected_coordinates_outer[i][1], node_coords[1], tolerance);
//         KRATOS_EXPECT_NEAR(expected_coordinates_outer[i][2], node_coords[2], tolerance);
//         i++;
//     }

// }



// // Tests the SnakeSbmProcess with an inner geometry
// KRATOS_TEST_CASE_IN_SUITE(SnakeSbmProcessNurbsInner, KratosIgaFastSuite)
// {
    
//     Model model;
//     ModelPart& iga_model_part = model.CreateModelPart("iga_model_part");

//     ModelPart& surrogate_sub_model_part_inner = iga_model_part.CreateSubModelPart("surrogate_inner");
//     iga_model_part.CreateSubModelPart("surrogate_outer");

//     ModelPart& skin_model_part_inner_initial = model.CreateModelPart("skin_model_part_inner_initial");
//     model.CreateModelPart("skin_model_part_outer_initial");

//     skin_model_part_inner_initial.CreateNewProperties(0);

//     // First arc of the circle
//     {
//         PointerVector<Node> control_points;
//         std::string condition_name = "SupportSolidCondition";
//         std::string layer_name = "left";

//         control_points.push_back(Node::Pointer(new Node(1, 0.0,
//                     0.9999999999999984,
//                     0.0)));
//         control_points.push_back(Node::Pointer(new Node(2, -1.7320508075688767,
//                     0.9999999999999984,
//                     0.0)));
//         control_points.push_back(Node::Pointer(new Node(3,  -0.8660254037844387,
//                     -0.4999999999999998,
//                     0.0)));
//         control_points.push_back(Node::Pointer(new Node(4,  -2.449293598294706e-16,
//                     -1.9999999999999996,
//                     0.0)));
//         control_points.push_back(Node::Pointer(new Node(5,  0.8660254037844384,
//                     -0.5000000000000004,
//                     0.0)));
//         control_points.push_back(Node::Pointer(new Node(6,  1.7320508075688776,
//                     0.9999999999999984,
//                     0.0)));
//         control_points.push_back(Node::Pointer(new Node(7,  2.4492935982947064e-16,
//                     0.9999999999999984,
//                     0.0)));

//         std::vector<double> weights_vector_temp = {1.0,
//                 0.5000000000000001,
//                 1.0,
//                 0.5000000000000001,
//                 1.0,
//                 0.5000000000000001,
//                 1.0};
//         int polynomial_degree = 2;
//         std::vector<double> knot_vector_temp = {0.0,
//                 0.0,
//                 0.0,
//                 0.3333333333333333,
//                 0.3333333333333333,
//                 0.6666666666666666,
//                 0.6666666666666666,
//                 1.0,
//                 1.0,
//                 1.0};

//         Vector weights_vector(weights_vector_temp.size());
//         for (std::size_t i = 0; i < weights_vector_temp.size(); ++i) 
//             weights_vector[i] = weights_vector_temp[i];
        
//         Vector knot_vector(knot_vector_temp.size());
//         for (std::size_t i = 0; i < knot_vector_temp.size(); ++i) 
//             knot_vector[i] = knot_vector_temp[i];

//         // Create the NURBS curve geometry
//         using NurbsCurveGeometryPointerType = NurbsCurveGeometry<2, PointerVector<Node>>::Pointer;
//         NurbsCurveGeometryPointerType p_curve(new NurbsCurveGeometry<2, PointerVector<Node>>(
//                                                             control_points,
//                                                             polynomial_degree,
//                                                             knot_vector, 
//                                                             weights_vector)); 
            
//         // link the boundary condition and layer name to the nurbs curve 
//         p_curve->SetValue(CONDITION_NAME, condition_name);
//         p_curve->SetValue(IDENTIFIER, layer_name);

//         p_curve->SetId(0);
//         skin_model_part_inner_initial.AddGeometry(p_curve);
//     }

//     ModelPart& skin_model_part = model.CreateModelPart("skin_model_part");
//     skin_model_part.CreateSubModelPart("inner");
//     skin_model_part.CreateSubModelPart("outer");
    
//     const std::vector<double> list_knot_u = {-2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0};
//     const std::vector<double> list_knot_v = {-2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0};

//     Vector unique_knot_vector_u(list_knot_u.size());
//     Vector unique_knot_vector_v(list_knot_v.size());

//     // Fill the vectors using a loop
//     for (std::size_t i = 0; i < list_knot_u.size(); ++i) 
//         unique_knot_vector_u[i] = list_knot_u[i];
//     for (std::size_t i = 0; i < list_knot_v.size(); ++i) 
//         unique_knot_vector_v[i] = list_knot_v[i];

//     Kratos::Parameters snake_parameters(R"(
//         {
//             "model_part_name" : "iga_model_part",
//             "skin_model_part_inner_initial_name" : "skin_model_part_inner_initial",
//             "skin_model_part_outer_initial_name" : "skin_model_part_outer_initial",
//             "skin_model_part_name" : "skin_model_part",
//             "echo_level" : 0,
//             "lambda_inner" : 0.5,
//             "number_of_inner_loops": 1
//         }
//     )");

//     SnakeSbmProcess snake_sbm_process(model, snake_parameters);

//     iga_model_part.SetValue(KNOT_VECTOR_U, unique_knot_vector_u);
//     iga_model_part.SetValue(KNOT_VECTOR_V, unique_knot_vector_v);

//     snake_sbm_process.Execute() ;
    
//     const double tolerance = 1.0e-6;

//     // Expected coordinates of nodes (modify according to actual expected values)
//     std::vector<std::array<double, 3>> expected_coordinates = {
//             {-0.8, -0.4, 0},
//             {0.8, -0.4, 0},
//             {-0.8, 0, 0},
//             {0.8, 0, 0},
//             {-0.8, 0.4, 0},
//             {0.8, 0.4, 0},
//             {-0.8, 0.8, 0},
//             {0.8, 0.8, 0},
//             {-0.8, -0.8, 0},
//             {-0.8, 0.8, 0},
//             {-0.4, -0.8, 0},
//             {-0.4, 0.8, 0},
//             {0, -0.8, 0},
//             {0, 0.8, 0},
//             {0.4, -0.8, 0},
//             {0.4, 0.8, 0}
//     };
    
//     // Ensure the number of nodes matches expectation
//     KRATOS_EXPECT_NEAR(surrogate_sub_model_part_inner.NumberOfConditions(), expected_coordinates.size(), tolerance);
    
//     // Iterate over nodes and compare coordinates
//     unsigned int i = 0;
//     for (auto& cond : surrogate_sub_model_part_inner.Conditions()) {
//         std::array<double, 3> node_coords = {cond.GetGeometry()[0].X(), cond.GetGeometry()[0].Y(), cond.GetGeometry()[0].Z()};

//         // Compare each coordinate
//         KRATOS_EXPECT_NEAR(expected_coordinates[i][0], node_coords[0], tolerance);
//         KRATOS_EXPECT_NEAR(expected_coordinates[i][1], node_coords[1], tolerance);
//         KRATOS_EXPECT_NEAR(expected_coordinates[i][2], node_coords[2], tolerance);
//         i++;
//     }
// }


}
