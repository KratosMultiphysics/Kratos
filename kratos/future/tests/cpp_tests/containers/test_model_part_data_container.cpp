//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//  Collaborator:    Vicente Mataix Ferrandiz
//
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "future/containers/model_part_data_container.h"
#include "geometries/triangle_2d_3.h"
#include "includes/expect.h"
#include "includes/variables.h"
#include "testing/testing.h"
#include "utilities/variable_utils.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(FutureModelPartDataContainerNodes, KratosCoreFutureSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("test");
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_model_part.SetBufferSize(3);
    auto p_node1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node2 = r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto p_node3 = r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);

    Future::ModelPartDataContainer model_part_data(r_model_part);
    auto& r_nodes_data = model_part_data.Nodes();

    // The snapshot registered every node with the model part's buffer size
    KRATOS_EXPECT_EQ(r_nodes_data.NumberOfEntities(), 3);
    KRATOS_EXPECT_EQ(r_nodes_data.GetBufferSize(), 3);
    KRATOS_EXPECT_TRUE(r_nodes_data.HasEntity(2));

    // Historical (TimeStep-buffered) and non-historical variables through the new path
    r_nodes_data.AddHistoricalVariable(TEMPERATURE);
    r_nodes_data.AddVariable(PRESSURE);

    for (auto& r_node : r_model_part.Nodes()) {
        r_nodes_data.SetValue(r_node, TEMPERATURE, 100.0 + r_node.Id());
        r_nodes_data.SetValue(r_node, PRESSURE, 200.0 + r_node.Id());
    }

    // Advance the buffer in lockstep with the legacy workflow
    r_model_part.CloneTimeStep(1.0);
    model_part_data.CloneStepData(StepCategory::TimeStep);
    r_nodes_data.SetValue(*p_node1, TEMPERATURE, 111.0);

    KRATOS_EXPECT_EQ(r_nodes_data.GetValue(*p_node1, TEMPERATURE), 111.0);
    KRATOS_EXPECT_EQ(r_nodes_data.GetValue(*p_node1, TEMPERATURE, 1), 101.0); // previous step
    KRATOS_EXPECT_EQ(r_nodes_data.GetValue(*p_node2, TEMPERATURE), 102.0);    // cloned
    KRATOS_EXPECT_EQ(r_nodes_data.GetValue(*p_node3, PRESSURE), 203.0);
}

KRATOS_TEST_CASE_IN_SUITE(FutureModelPartDataContainerNoCrossContamination, KratosCoreFutureSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("test");
    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    r_model_part.SetBufferSize(2);
    auto p_node = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);

    // Seed the LEGACY storage first
    p_node->FastGetSolutionStepValue(TEMPERATURE) = 25.0; // historical
    p_node->SetValue(PRESSURE, 5.0);                      // non-historical

    Future::ModelPartDataContainer model_part_data(r_model_part);
    auto& r_nodes_data = model_part_data.Nodes();
    r_nodes_data.AddHistoricalVariable(TEMPERATURE);
    r_nodes_data.AddVariable(PRESSURE);

    // The new path starts at zero regardless of the legacy values
    KRATOS_EXPECT_EQ(r_nodes_data.GetValue(*p_node, TEMPERATURE), 0.0);
    KRATOS_EXPECT_EQ(r_nodes_data.GetValue(*p_node, PRESSURE), 0.0);

    // Writing through the new path leaves the legacy storage untouched...
    r_nodes_data.SetValue(*p_node, TEMPERATURE, 999.0);
    r_nodes_data.SetValue(*p_node, PRESSURE, 888.0);
    KRATOS_EXPECT_EQ(p_node->FastGetSolutionStepValue(TEMPERATURE), 25.0);
    KRATOS_EXPECT_EQ(p_node->GetValue(PRESSURE), 5.0);

    // ...and writing through the legacy storage leaves the new path untouched
    p_node->FastGetSolutionStepValue(TEMPERATURE) = 26.0;
    p_node->SetValue(PRESSURE, 6.0);
    KRATOS_EXPECT_EQ(r_nodes_data.GetValue(*p_node, TEMPERATURE), 999.0);
    KRATOS_EXPECT_EQ(r_nodes_data.GetValue(*p_node, PRESSURE), 888.0);
}

KRATOS_TEST_CASE_IN_SUITE(FutureModelPartDataContainerUpdate, KratosCoreFutureSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("test");
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);

    Future::ModelPartDataContainer model_part_data(r_model_part);
    model_part_data.Nodes().AddVariable(PRESSURE);
    model_part_data.Nodes().SetValue(std::size_t(1), PRESSURE, 1.5);
    KRATOS_EXPECT_EQ(model_part_data.Nodes().NumberOfEntities(), 1);

    // Nodes created after the snapshot are unknown until Update()
    auto p_new_node = r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    KRATOS_EXPECT_FALSE(model_part_data.Nodes().HasEntity(2));

    model_part_data.Update();

    KRATOS_EXPECT_TRUE(model_part_data.Nodes().HasEntity(2));
    KRATOS_EXPECT_EQ(model_part_data.Nodes().GetValue(*p_new_node, PRESSURE), 0.0);
    KRATOS_EXPECT_EQ(model_part_data.Nodes().GetValue(std::size_t(1), PRESSURE), 1.5); // preserved across growth
}

KRATOS_TEST_CASE_IN_SUITE(FutureModelPartDataContainerElementsAndConditions, KratosCoreFutureSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("test");
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    auto p_prop = r_model_part.CreateNewProperties(1);

    // Two elements SHARING the same three nodes, plus one condition
    auto p_elem1 = r_model_part.CreateNewElement("Element2D3N", 1, {{1, 2, 3}}, p_prop);
    auto p_elem2 = r_model_part.CreateNewElement("Element2D3N", 2, {{1, 2, 3}}, p_prop);
    auto p_cond = r_model_part.CreateNewCondition("SurfaceCondition3D3N", 1, {{1, 2, 3}}, p_prop);

    Future::ModelPartDataContainer model_part_data(r_model_part);
    auto& r_elements_data = model_part_data.Elements();
    auto& r_conditions_data = model_part_data.Conditions();

    KRATOS_EXPECT_EQ(r_elements_data.NumberOfEntities(), 2);
    KRATOS_EXPECT_EQ(r_conditions_data.NumberOfEntities(), 1);

    // New path: values are keyed by element Id, so the two elements are INDEPENDENT
    r_elements_data.AddVariable(PRESSURE);
    r_elements_data.SetValue(*p_elem1, PRESSURE, 1.0);
    r_elements_data.SetValue(*p_elem2, PRESSURE, 2.0);
    KRATOS_EXPECT_EQ(r_elements_data.GetValue(*p_elem1, PRESSURE), 1.0);
    KRATOS_EXPECT_EQ(r_elements_data.GetValue(*p_elem2, PRESSURE), 2.0);

    // Legacy path is untouched by the new-path writes and keeps its own semantics
    // (elements DO share their geometry data when they share a geometry object; these two
    // elements have distinct Triangle2D3 geometries over the same nodes, so SetValue is
    // per element here as well - the important check is the mutual independence)
    p_elem1->SetValue(PRESSURE, 10.0);
    KRATOS_EXPECT_EQ(p_elem1->GetValue(PRESSURE), 10.0);
    KRATOS_EXPECT_EQ(r_elements_data.GetValue(*p_elem1, PRESSURE), 1.0); // new path unaffected

    // Conditions behave the same way
    r_conditions_data.AddVariable(TEMPERATURE);
    r_conditions_data.SetValue(*p_cond, TEMPERATURE, 42.0);
    KRATOS_EXPECT_EQ(r_conditions_data.GetValue(*p_cond, TEMPERATURE), 42.0);
    KRATOS_EXPECT_FALSE(p_cond->Has(TEMPERATURE)); // legacy storage untouched
}

KRATOS_TEST_CASE_IN_SUITE(FutureModelPartDataContainerGeometries, KratosCoreFutureSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("test");
    auto p_node1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node2 = r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto p_node3 = r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    r_model_part.CreateNewGeometry("Triangle2D3", 1, {{1, 2, 3}});
    r_model_part.CreateNewGeometry("Triangle2D3", 2, {{1, 2, 3}});

    Future::ModelPartDataContainer model_part_data(r_model_part);
    auto& r_geometries_data = model_part_data.Geometries();

    KRATOS_EXPECT_EQ(r_geometries_data.NumberOfEntities(), 2);

    r_geometries_data.AddVariable(DISTANCE);
    auto& r_geom1 = r_model_part.GetGeometry(1);
    r_geometries_data.SetValue(r_geom1, DISTANCE, 3.5);
    KRATOS_EXPECT_EQ(r_geometries_data.GetValue(r_geom1, DISTANCE), 3.5);
    KRATOS_EXPECT_EQ(r_geometries_data.GetValue(std::size_t(2), DISTANCE), 0.0);

    // The legacy geometry data container is untouched
    KRATOS_EXPECT_FALSE(r_geom1.Has(DISTANCE));

    // Standalone geometries (not in a ModelPart) work against a raw EntityDataContainer,
    // provided they carry a meaningful unique Id
    auto standalone_geometry = Triangle2D3<Node>(p_node1, p_node2, p_node3);
    standalone_geometry.SetId(77);

    Future::EntityDataContainer standalone_data(1, 2);
    standalone_data.RegisterEntity(standalone_geometry);
    standalone_data.AddVariable(TEMPERATURE);
    standalone_data.SetValue(standalone_geometry, TEMPERATURE, 7.5);
    KRATOS_EXPECT_EQ(standalone_data.GetValue(standalone_geometry, TEMPERATURE), 7.5);
}

KRATOS_TEST_CASE_IN_SUITE(FutureModelPartDataContainerMasterSlaveConstraints, KratosCoreFutureSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("test");
    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(REACTION);
    auto p_node1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node2 = r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X, REACTION_X);
    }
    auto p_constraint = r_model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 1, *p_node1, DISPLACEMENT_X, *p_node2, DISPLACEMENT_X, 1.0, 0);

    Future::ModelPartDataContainer model_part_data(r_model_part);
    auto& r_constraints_data = model_part_data.MasterSlaveConstraints();

    KRATOS_EXPECT_EQ(r_constraints_data.NumberOfEntities(), 1);
    KRATOS_EXPECT_TRUE(r_constraints_data.HasEntity(1));

    // Read/write through the new path
    r_constraints_data.AddVariable(DISTANCE);
    r_constraints_data.SetValue(*p_constraint, DISTANCE, 12.5);
    KRATOS_EXPECT_EQ(r_constraints_data.GetValue(*p_constraint, DISTANCE), 12.5);

    // No cross-contamination with the constraint's legacy DataValueContainer
    KRATOS_EXPECT_FALSE(p_constraint->Has(DISTANCE));
    p_constraint->SetValue(DISTANCE, 99.0);
    KRATOS_EXPECT_EQ(p_constraint->GetValue(DISTANCE), 99.0);
    KRATOS_EXPECT_EQ(r_constraints_data.GetValue(*p_constraint, DISTANCE), 12.5);
}

KRATOS_TEST_CASE_IN_SUITE(FutureModelPartDataContainerEndToEnd, KratosCoreFutureSuite)
{
    // Realistic usage: add a new-path variable, loop over the nodes writing through the
    // bulk span, read back per entity
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("test");
    for (std::size_t id = 1; id <= 10; ++id) {
        r_model_part.CreateNewNode(id, 1.0 * id, 0.0, 0.0);
    }

    Future::ModelPartDataContainer model_part_data(r_model_part);
    auto& r_nodes_data = model_part_data.Nodes();
    auto accessor = r_nodes_data.AddVariable(DISTANCE);

    auto distance_span = r_nodes_data.GetDataSpan(accessor);
    for (const auto& r_node : r_model_part.Nodes()) {
        distance_span[r_nodes_data.Index(r_node)] = r_node.X();
    }

    for (const auto& r_node : r_model_part.Nodes()) {
        KRATOS_EXPECT_EQ(r_nodes_data.GetValue(r_node, DISTANCE), r_node.X());
    }
}

} // namespace Kratos::Testing
