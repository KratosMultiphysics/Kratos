//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// Project includes
#include "includes/checks.h"
#include "testing/testing.h"
#include "modeler/create_entities_from_geometries_modeler.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CreateEntitiesFromGeometriesModeler, KratosCoreFastSuite)
{
    // Set up the test model
    Model model;
    auto& r_model_part_tri = model.CreateModelPart("Triangles");
    auto& r_model_part_tri_2d = r_model_part_tri.CreateSubModelPart("Tri2D");
    auto& r_model_part_tri_3d = r_model_part_tri.CreateSubModelPart("Tri3D");
    auto& r_model_part_quads = model.CreateModelPart("Quadrilaterals");
    auto& r_model_part_quads_1 = r_model_part_quads.CreateSubModelPart("Quads1");
    auto& r_model_part_quads_2 = r_model_part_quads.CreateSubModelPart("Quads2");
    auto& r_model_part_quads_boundary = r_model_part_quads.CreateSubModelPart("QuadsBoundary");

    r_model_part_tri.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part_tri.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_model_part_tri.CreateNewNode(3, 0.0, 1.0, 0.0);
    r_model_part_tri_2d.CreateNewGeometry("Triangle2D3", 1, {{1, 2, 3}});
    r_model_part_tri_3d.CreateNewGeometry("Triangle3D3", 2, {{1, 2, 3}});

    r_model_part_quads.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part_quads.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_model_part_quads.CreateNewNode(3, 1.0, 1.0, 0.0);
    r_model_part_quads.CreateNewNode(4, 0.0, 1.0, 0.0);
    r_model_part_quads.CreateNewNode(5, 2.0, 0.0, 0.0);
    r_model_part_quads.CreateNewNode(6, 2.0, 1.0, 0.0);
    r_model_part_quads_1.CreateNewGeometry("Quadrilateral2D4", 1, {{1, 2, 3, 4}});
    r_model_part_quads_2.CreateNewGeometry("Quadrilateral2D4", 2, {{2, 5, 6, 3}});
    r_model_part_quads_boundary.CreateNewGeometry("Line2D2", 3, {{1, 2}});
    r_model_part_quads_boundary.CreateNewGeometry("Line2D2", 4, {{5, 6}});

    // Create and execute the modeler
    Parameters settings(R"({
        "elements_list" : [{
            "model_part_name" : "Triangles.Tri2D",
            "element_name" : "Element2D3N"
        },{
            "model_part_name" : "Quadrilaterals.Quads1",
            "element_name" : "Element2D4N"
        },{
            "model_part_name" : "Quadrilaterals.Quads2",
            "element_name" : "Element2D4N"
        }],
        "conditions_list" : [{
            "model_part_name" : "Quadrilaterals.QuadsBoundary",
            "condition_name" : "LineCondition2D2N"
        },{
            "model_part_name" : "Triangles.Tri3D",
            "condition_name" : "SurfaceCondition3D3N"
        }]
    })");
    auto modeler = CreateEntitiesFromGeometriesModeler(model, settings);
    modeler.SetupGeometryModel();
    modeler.PrepareGeometryModel();
    modeler.SetupModelPart();

    // Check results
    KRATOS_CHECK_EQUAL(r_model_part_tri.NumberOfElements(), 1);
    KRATOS_CHECK_EQUAL(r_model_part_tri.NumberOfConditions(), 1);
    KRATOS_CHECK_EQUAL(r_model_part_tri_2d.NumberOfElements(), 1);
    KRATOS_CHECK_EQUAL(r_model_part_tri_2d.NumberOfConditions(), 0);
    KRATOS_CHECK_EQUAL(r_model_part_tri_3d.NumberOfElements(), 0);
    KRATOS_CHECK_EQUAL(r_model_part_tri_3d.NumberOfConditions(), 1);
    KRATOS_CHECK_EQUAL(r_model_part_quads.NumberOfElements(), 2);
    KRATOS_CHECK_EQUAL(r_model_part_quads.NumberOfConditions(), 2);
    KRATOS_CHECK_EQUAL(r_model_part_quads_1.NumberOfElements(), 1);
    KRATOS_CHECK_EQUAL(r_model_part_quads_1.NumberOfConditions(), 0);
    KRATOS_CHECK_EQUAL(r_model_part_quads_2.NumberOfElements(), 1);
    KRATOS_CHECK_EQUAL(r_model_part_quads_2.NumberOfConditions(), 0);
    KRATOS_CHECK_EQUAL(r_model_part_quads_boundary.NumberOfElements(), 0);
    KRATOS_CHECK_EQUAL(r_model_part_quads_boundary.NumberOfConditions(), 2);
}

KRATOS_TEST_CASE_IN_SUITE(CreateEntitiesFromGeometriesModelerMixedEntities, KratosCoreFastSuite)
{
    // Set up the test model
    Model model;
    auto& r_model_part = model.CreateModelPart("MainModelPart");

    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
    r_model_part.CreateNewNode(4, 0.0, 1.0, 0.0);
    r_model_part.CreateNewNode(5, 2.0, 0.0, 0.0);
    r_model_part.CreateNewNode(6, 2.0, 1.0, 0.0);
    r_model_part.CreateNewGeometry("Line2D2", 1, {{1, 2}});
    r_model_part.CreateNewGeometry("Line2D2", 2, {{5, 6}});
    r_model_part.CreateNewGeometry("Triangle2D3", 3, {1, 2, 4});
    r_model_part.CreateNewGeometry("Quadrilateral2D4", 4, {1, 2, 3, 4});
    r_model_part.CreateNewGeometry("Quadrilateral2D4", 5, {2, 5, 6, 3});

    // Create and execute the modeler
    Parameters settings(R"({
        "elements_list" : [{
            "model_part_name" : "MainModelPart",
            "element_name" : "Element2D3N;Element2D4N"
        }],
        "conditions_list" : [{
            "model_part_name" : "MainModelPart",
            "condition_name" : "LineCondition2D2N"
        }]
    })");
    auto modeler = CreateEntitiesFromGeometriesModeler(model, settings);
    modeler.SetupGeometryModel();
    modeler.PrepareGeometryModel();
    modeler.SetupModelPart();

    // Check results
    KRATOS_CHECK_EQUAL(r_model_part.NumberOfElements(), 3);
    KRATOS_CHECK_EQUAL(r_model_part.NumberOfConditions(), 2);
}

}  // namespace Kratos::Testing.
