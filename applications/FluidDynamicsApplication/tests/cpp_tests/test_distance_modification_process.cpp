//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "tests/cpp_tests/fluid_dynamics_fast_suite.h"
#include "custom_processes/distance_modification_process.h"

namespace Kratos {
namespace Testing {

void TriangleModelPartForDistanceModification(
    const bool ContinuousDistance,
    ModelPart& rModelPart) {

    rModelPart.SetBufferSize(3);
    rModelPart.AddNodalSolutionStepVariable(DISTANCE);
    Properties::Pointer p_properties = rModelPart.CreateNewProperties(0);

    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> element_nodes{1, 2, 3};
    rModelPart.CreateNewElement("Element2D3N", 1, element_nodes, p_properties);

    // Distance data
    if (ContinuousDistance){
        auto& r_geometry = (rModelPart.ElementsBegin())->GetGeometry();
        r_geometry[0].FastGetSolutionStepValue(DISTANCE) = -1.0e-5;
        r_geometry[1].FastGetSolutionStepValue(DISTANCE) = 1.0;
        r_geometry[2].FastGetSolutionStepValue(DISTANCE) = 1.0;
    } else {
        Vector elem_dist(3,1.0);
        elem_dist(0) = -1e-5;
        rModelPart.ElementsBegin()->SetValue(ELEMENTAL_DISTANCES, elem_dist);
        Vector elem_edge_dist(3,-1.0);
        elem_edge_dist(1) = 0.8;
        rModelPart.ElementsBegin()->SetValue(ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED, elem_edge_dist);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DistanceModificationTriangle, FluidDynamicsApplicationFastSuite) {
    Model model;
    ModelPart& model_part = model.CreateModelPart("TestPart");
    TriangleModelPartForDistanceModification(true, model_part);

    Parameters distance_mod_params( R"(
    {
        "model_part_name"                             : "TestPart",
        "distance_threshold"                          : 0.001,
        "continuous_distance"                         : true,
        "avoid_almost_empty_elements"                 : true,
        "deactivate_full_negative_elements"           : true,
        "recover_original_distance_at_each_step"      : true,
        "full_negative_elements_fixed_variables_list" : []
    }  )" );

    DistanceModificationProcess dist_mod_process(model_part,distance_mod_params);
    // Check the distance field modification
    dist_mod_process.ExecuteInitialize();
    dist_mod_process.ExecuteBeforeSolutionLoop();
    dist_mod_process.ExecuteInitializeSolutionStep();
    const double tolerance = 1e-9;
    const std::array<double, 3> expected_values = {-1.0e-3, 1.0, 1.0};
    for (size_t i_node = 1; i_node < 4; ++i_node) {
        KRATOS_EXPECT_NEAR((model_part.GetNode(i_node)).FastGetSolutionStepValue(DISTANCE), expected_values[i_node - 1], tolerance);
    }

    // Check that the flag TO_SPLIT is correctly set
    KRATOS_EXPECT_TRUE((model_part.ElementsBegin())->Is(TO_SPLIT));

    // Check the original distance recovering
    dist_mod_process.ExecuteFinalizeSolutionStep();
    const std::array<double, 3> expected_orig_values = {-1.0e-5, 1.0, 1.0};
    for (size_t i_node = 1; i_node < 4; ++i_node) {
        KRATOS_EXPECT_NEAR((model_part.GetNode(i_node)).FastGetSolutionStepValue(DISTANCE), expected_orig_values[i_node - 1], tolerance);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceModificationTriangle, FluidDynamicsApplicationFastSuite) {
    Model model;
    ModelPart& model_part = model.CreateModelPart("TestPart");
    TriangleModelPartForDistanceModification(false, model_part);

    Parameters distance_mod_params( R"(
    {
        "model_part_name"                        : "TestPart",
        "distance_threshold"                     : 0.001,
        "continuous_distance"                    : false,
        "avoid_almost_empty_elements"            : false,
        "deactivate_full_negative_elements"      : false,
        "recover_original_distance_at_each_step" : true
    }  )" );

    DistanceModificationProcess dist_mod_process(model_part,distance_mod_params);
    // Check the distance field modification
    dist_mod_process.ExecuteInitialize();
    dist_mod_process.ExecuteBeforeSolutionLoop();
    dist_mod_process.ExecuteInitializeSolutionStep();
    auto& r_elem_dist = (model_part.ElementsBegin())->GetValue(ELEMENTAL_DISTANCES);
    auto& r_elem_edge_dist_extra = (model_part.ElementsBegin())->GetValue(ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED);
    const double tolerance = 1e-9;
    const std::vector<double> expected_values = {-0.000797885, 1.0, 1.0};
    const std::vector<double> expected_extra_values = {-1.0, 0.9992027516, -1.0};
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist, expected_values, tolerance);
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_edge_dist_extra, expected_extra_values, tolerance);

    // Check that the flag TO_SPLIT is correctly set
    KRATOS_EXPECT_TRUE((model_part.ElementsBegin())->Is(TO_SPLIT));

    // Check the original distance recovering
    dist_mod_process.ExecuteFinalizeSolutionStep();
    r_elem_dist = (model_part.ElementsBegin())->GetValue(ELEMENTAL_DISTANCES);
    r_elem_edge_dist_extra = (model_part.ElementsBegin())->GetValue(ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED);
    const std::vector<double> expected_orig_values = {-1.0e-5, 1.0, 1.0};
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist, expected_orig_values, tolerance);
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_edge_dist_extra, expected_extra_values, tolerance);
}

}
}
