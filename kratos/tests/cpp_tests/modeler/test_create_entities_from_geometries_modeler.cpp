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
#include "testing/testing.h"
#include "modeler/create_entities_from_geometries_modeler.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CreateEntitiesFromGeometriesModeler, KratosCoreFastSuite)
{
    // Set up the test model
    Model model;
    auto& r_model_part_tri = model.CreateModelPart("Triangles");
    auto& r_model_part_quads = model.CreateModelPart("Quadrilaterals");
    auto& r_quads_1 = r_model_part_quads.CreateSubModelPart("Quads1");
    auto& r_quads_2 = r_model_part_quads.CreateSubModelPart("Quads2");
    auto& r_quads_boundary = r_model_part_quads.CreateSubModelPart("QuadsBoundary");

    r_model_part_tri.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part_tri.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_model_part_tri.CreateNewNode(3, 0.0, 1.0, 0.0);
    r_model_part_tri.CreateNewGeometry("Triangle2D3", 1, {1, 2, 3});

    r_model_part_quads.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part_quads.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_model_part_quads.CreateNewNode(3, 1.0, 1.0, 0.0);
    r_model_part_quads.CreateNewNode(4, 0.0, 1.0, 0.0);
    r_model_part_quads.CreateNewNode(5, 2.0, 0.0, 0.0);
    r_model_part_quads.CreateNewNode(6, 2.0, 1.0, 0.0);
    r_model_part_quads.CreateNewGeometry("Quadrilateral2D4", 1, {1, 2, 3, 4});
    r_model_part_quads.CreateNewGeometry("Quadrilateral2D4", 2, {2, 5, 6, 3});
    r_model_part_quads.CreateNewGeometry("Line2D2", 3, {{1, 2}});
    r_model_part_quads.CreateNewGeometry("Line2D2", 4, {{5, 6}});

    // Create the modeler
    Parameters settings(R"({
        "elements_list" : [{
            "model_part_name" : "Triangles",
            "element_name" : "Element2D3N"
        }],
        "conditions_list" : []
    })");
    auto modeler = CreateEntitiesFromGeometriesModeler(model, settings);
    modeler.SetupGeometryModel();
    modeler.PrepareGeometryModel();
    modeler.SetupModelPart();

    // KRATOS_CHECK_EQUAL(model_part_1.NumberOfProperties(), 2);
    // KRATOS_CHECK_EQUAL(model_part_2.NumberOfProperties(), 2);
    // KRATOS_CHECK_EQUAL(model_part_2.GetElement(1).GetProperties().Id(), 1);

    // model_part_2.GetProperties(1).SetValue(DISTANCE, 2.2);
    // KRATOS_CHECK_EQUAL(model_part_1.GetProperties(1).GetValue(DISTANCE), 1.1);
    // KRATOS_CHECK_EQUAL(model_part_2.GetProperties(1).GetValue(DISTANCE), 2.2);
}

}  // namespace Kratos::Testing.
