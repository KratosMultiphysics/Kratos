//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//                   Vicente Mataix Ferrandiz
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "spatial_containers/geometrical_objects_bins.h"
#include "containers/model.h"
#include "geometries/triangle_3d_3.h"

namespace Kratos::Testing 
{

ModelPart& CreateCubeSkinModelPart(
    Model& rCurrentModel,
    const double HalfX = 0.6,
    const double HalfY = 0.9,
    const double HalfZ = 0.3
    )
{
    // Generate the cube skin
    ModelPart& r_skin_part = rCurrentModel.CreateModelPart("Skin");
    r_skin_part.CreateNewNode(1, -HalfX, -HalfY, -HalfZ);
    r_skin_part.CreateNewNode(2,  HalfX, -HalfY, -HalfZ);
    r_skin_part.CreateNewNode(3,  HalfX,  HalfY, -HalfZ);
    r_skin_part.CreateNewNode(4, -HalfX,  HalfY, -HalfZ);
    r_skin_part.CreateNewNode(5, -HalfX, -HalfY,  HalfZ);
    r_skin_part.CreateNewNode(6,  HalfX, -HalfY,  HalfZ);
    r_skin_part.CreateNewNode(7,  HalfX,  HalfY,  HalfZ);
    r_skin_part.CreateNewNode(8, -HalfX,  HalfY,  HalfZ);
    Properties::Pointer p_properties(new Properties(0));
    r_skin_part.CreateNewElement("Element3D3N",  1, { 1,2,3 }, p_properties);
    r_skin_part.CreateNewElement("Element3D3N",  2, { 1,3,4 }, p_properties);
    r_skin_part.CreateNewElement("Element3D3N",  3, { 5,6,7 }, p_properties);
    r_skin_part.CreateNewElement("Element3D3N",  4, { 5,7,8 }, p_properties);
    r_skin_part.CreateNewElement("Element3D3N",  5, { 3,6,2 }, p_properties);
    r_skin_part.CreateNewElement("Element3D3N",  6, { 3,7,6 }, p_properties);
    r_skin_part.CreateNewElement("Element3D3N",  7, { 4,5,1 }, p_properties);
    r_skin_part.CreateNewElement("Element3D3N",  8, { 4,8,5 }, p_properties);
    r_skin_part.CreateNewElement("Element3D3N",  9, { 3,4,8 }, p_properties);
    r_skin_part.CreateNewElement("Element3D3N", 10, { 3,8,7 }, p_properties);
    r_skin_part.CreateNewElement("Element3D3N", 11, { 2,1,5 }, p_properties);
    r_skin_part.CreateNewElement("Element3D3N", 12, { 2,5,6 }, p_properties);

    return r_skin_part;
}

ModelPart& CreateCubeModelPart(Model& rCurrentModel)
{
    // Generate the cube skin
    ModelPart& r_model_part = rCurrentModel.CreateModelPart("Cube");

    // Create properties
    auto p_properties = r_model_part.CreateNewProperties(1, 0);

    r_model_part.CreateNewNode(1 , 0.0 , 1.0 , 1.0);
    r_model_part.CreateNewNode(2 , 0.0 , 1.0 , 0.0);
    r_model_part.CreateNewNode(3 , 0.0 , 0.0 , 1.0);
    r_model_part.CreateNewNode(4 , 0.5 , 1.0 , 1.0);
    r_model_part.CreateNewNode(5 , 0.0 , 0.0 , 0.0);
    r_model_part.CreateNewNode(6 , 0.5 , 1.0 , 0.0);
    r_model_part.CreateNewNode(7 , 0.5 , 0.0 , 1.0);
    r_model_part.CreateNewNode(8 , 0.5 , 0.0 , 0.0);
    r_model_part.CreateNewNode(9 , 1.0 , 1.0 , 1.0);
    r_model_part.CreateNewNode(10 , 1.0 , 1.0 , 0.0);
    r_model_part.CreateNewNode(11 , 1.0 , 0.0 , 1.0);
    r_model_part.CreateNewNode(12 , 1.0 , 0.0 , 0.0);

    // Create elements
    r_model_part.CreateNewElement("Element3D8N", 1, {{5,8,6,2,3,7,4,1}}, p_properties);
    r_model_part.CreateNewElement("Element3D8N", 2, {{8,12,10,6,7,11,9,4}}, p_properties);

    return r_model_part;
}

/** Checks bins bounding box
*/
KRATOS_TEST_CASE_IN_SUITE(GeometricalObjectsBinsBoundingBox, KratosCoreFastSuite) 
{
    constexpr double tolerance = 1e-12;

    Model current_model;

    const double cube_x = 0.6;
    const double cube_y = 0.9;
    const double cube_z = 0.3;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeSkinModelPart(current_model, cube_x, cube_y, cube_z);

    GeometricalObjectsBins bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd());

    auto bounding_box = bins.GetBoundingBox();

    KRATOS_CHECK_NEAR(bounding_box.GetMinPoint()[0],-cube_x, tolerance);
    KRATOS_CHECK_NEAR(bounding_box.GetMinPoint()[1],-cube_y, tolerance);
    KRATOS_CHECK_NEAR(bounding_box.GetMinPoint()[2],-cube_z, tolerance);
    KRATOS_CHECK_NEAR(bounding_box.GetMaxPoint()[0], cube_x, tolerance);
    KRATOS_CHECK_NEAR(bounding_box.GetMaxPoint()[1], cube_y, tolerance);
    KRATOS_CHECK_NEAR(bounding_box.GetMaxPoint()[2], cube_z, tolerance);
}

/** Checks bins number of cells
*/
KRATOS_TEST_CASE_IN_SUITE(GeometricalObjectsBinsCellSizes, KratosCoreFastSuite) 
{
    constexpr double tolerance = 1e-12;

    Model current_model;

    const double cube_x = 0.6;
    const double cube_y = 0.9;
    const double cube_z = 0.3;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeSkinModelPart(current_model, cube_x, cube_y, cube_z);

    GeometricalObjectsBins bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd());

    auto number_of_cells = bins.GetNumberOfCells();
    KRATOS_CHECK_EQUAL(number_of_cells[0], 3);
    KRATOS_CHECK_EQUAL(number_of_cells[1], 3);
    KRATOS_CHECK_EQUAL(number_of_cells[2], 2);

    auto cell_sizes = bins.GetCellSizes();
    KRATOS_CHECK_NEAR(cell_sizes[0], 2.00 * cube_x / 3.00, tolerance);
    KRATOS_CHECK_NEAR(cell_sizes[1], 2.00 * cube_y / 3.00, tolerance);
    KRATOS_CHECK_NEAR(cell_sizes[2], 2.00 * cube_z / 2.00, tolerance);
}


/** Checks bins AddObjectsToCells
*/
KRATOS_TEST_CASE_IN_SUITE(GeometricalObjectsBinsAddObjectsToCells, KratosCoreFastSuite) 
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeSkinModelPart(current_model);

    GeometricalObjectsBins bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd());

    auto& cell = bins.GetCell(0,0,0);
    KRATOS_CHECK_EQUAL(cell.size(), 4);
    for(auto& geometrical_object : cell){
        std::size_t id = geometrical_object->GetId();
        KRATOS_CHECK((id == 1) ||(id == 2) ||(id == 7) ||(id == 11));
    }

    cell = bins.GetCell(2,2,1);
    KRATOS_CHECK_EQUAL(cell.size(), 4);
    for(auto& geometrical_object : cell){
        std::size_t id = geometrical_object->GetId();
        KRATOS_CHECK((id == 3) ||(id == 4) ||(id == 6) ||(id == 10));
    }
}

/** Checks bins search in radius
*/
KRATOS_TEST_CASE_IN_SUITE(GeometricalObjectsBinsSearchInRadius, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeSkinModelPart(current_model);

    GeometricalObjectsBins bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd());

    std::vector<GeometricalObjectsBins::ResultType> results;
    Point center_point{0.00,0.00,0.00};

    bins.SearchInRadius(center_point, .29, results);
    KRATOS_CHECK_EQUAL(results.size(), 0);

    bins.SearchInRadius(center_point, .3, results);
    KRATOS_CHECK_EQUAL(results.size(), 4);

    bins.SearchInRadius(center_point, .4, results);
    KRATOS_CHECK_EQUAL(results.size(), 4);

    bins.SearchInRadius(center_point, .6, results);
    KRATOS_CHECK_EQUAL(results.size(), 8);

    bins.SearchInRadius(center_point, .7, results);
    KRATOS_CHECK_EQUAL(results.size(), 8);

    bins.SearchInRadius(center_point, .9, results);
    KRATOS_CHECK_EQUAL(results.size(), 12);
}

/** Checks bins search in radius
*/
KRATOS_TEST_CASE_IN_SUITE(GeometricalObjectsBinsSearchInRadiusContainer, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeSkinModelPart(current_model);

    GeometricalObjectsBins bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd());

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    auto p_node = r_point_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto& r_array_nodes = r_point_model_part.Nodes();

    GeometricalObjectsBins::ResultTypeContainerMap results;

    // 0.29 radius
    bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.29, results);
    KRATOS_CHECK_EQUAL(results.NumberOfSearchResults(), 1);
    KRATOS_CHECK_IS_FALSE(results[*p_node].IsObjectFound());
    KRATOS_CHECK_EQUAL(results[*p_node].NumberOfGlobalResults(), 0);

    // 0.3 radius
    bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.3, results);
    KRATOS_CHECK_EQUAL(results.NumberOfSearchResults(), 1);
    KRATOS_CHECK(results[*p_node].IsObjectFound());
    KRATOS_CHECK_EQUAL(results[*p_node].NumberOfGlobalResults(), 4);

    // 0.4 radius
    bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.4, results);
    KRATOS_CHECK_EQUAL(results.NumberOfSearchResults(), 1);
    KRATOS_CHECK(results[*p_node].IsObjectFound());
    KRATOS_CHECK_EQUAL(results[*p_node].NumberOfGlobalResults(), 4);

    // 0.6 radius
    bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.6, results);
    KRATOS_CHECK_EQUAL(results.NumberOfSearchResults(), 1);
    KRATOS_CHECK(results[*p_node].IsObjectFound());
    KRATOS_CHECK_EQUAL(results[*p_node].NumberOfGlobalResults(), 8);

    // 0.7 radius
    bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.7, results);
    KRATOS_CHECK_EQUAL(results.NumberOfSearchResults(), 1);
    KRATOS_CHECK(results[*p_node].IsObjectFound());
    KRATOS_CHECK_EQUAL(results[*p_node].NumberOfGlobalResults(), 8);

    // 0.9 radius
    bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.9, results);
    KRATOS_CHECK_EQUAL(results.NumberOfSearchResults(), 1);
    KRATOS_CHECK(results[*p_node].IsObjectFound());
    KRATOS_CHECK_EQUAL(results[*p_node].NumberOfGlobalResults(), 12);
}

/** Checks bins search nearest
*/
KRATOS_TEST_CASE_IN_SUITE(GeometricalObjectsBinsSearchNearestInRadius, KratosCoreFastSuite) 
{
    constexpr double tolerance = 1e-12;

    Model current_model;

    // Cube coordinates
    const double cube_z = 0.3;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeSkinModelPart(current_model, 0.6, 0.9, cube_z);

    GeometricalObjectsBins bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd());

    double epsilon = 1.00e-6;
    Point near_point{epsilon,epsilon,epsilon};
    auto result = bins.SearchNearestInRadius(near_point, cube_z - 1.e-4);

    KRATOS_CHECK_IS_FALSE(result.IsObjectFound());

    result = bins.SearchNearestInRadius(near_point, cube_z + 1.e-4);
    KRATOS_CHECK_NEAR(result.GetDistance(), (cube_z - epsilon), tolerance);

    std::size_t id = result.Get()->Id();
    KRATOS_CHECK(id == 3);
}

/** Checks bins search nearest
*/
KRATOS_TEST_CASE_IN_SUITE(GeometricalObjectsBinsSearchNearestInRadiusContainer, KratosCoreFastSuite) 
{
    constexpr double tolerance = 1e-12;

    Model current_model;

    // Cube coordinates
    const double cube_z = 0.3;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeSkinModelPart(current_model, 0.6, 0.9, cube_z);

    GeometricalObjectsBins bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd());

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    const double epsilon = 1.0e-6;
    auto p_node = r_point_model_part.CreateNewNode(1, epsilon, epsilon, epsilon);
    auto& r_array_nodes = r_point_model_part.Nodes();

    GeometricalObjectsBins::ResultTypeContainerMap results;
    bins.SearchNearestInRadius(r_array_nodes.begin(), r_array_nodes.end(), cube_z - 1.e-4, results);

    KRATOS_CHECK_EQUAL(results.NumberOfSearchResults(), 1);
    KRATOS_CHECK_IS_FALSE(results[*p_node].IsObjectFound());

    bins.SearchNearestInRadius(r_array_nodes.begin(), r_array_nodes.end(), cube_z + 1.e-4, results);

    KRATOS_CHECK_EQUAL(results.NumberOfSearchResults(), 1);
    KRATOS_CHECK(results[*p_node].IsObjectFound());
    KRATOS_CHECK_EQUAL(results[*p_node].NumberOfGlobalResults(), 1);

    // Distances
    auto distances = results[*p_node].GetDistances();
    KRATOS_CHECK_NEAR(distances[0], (cube_z - epsilon), tolerance);

    // Compute indices
    auto indices = results[*p_node].GetResultIndices();
    const std::size_t id = indices[0];
    KRATOS_CHECK_EQUAL(id, 3);
}

/** Checks bins search nearest
*/
KRATOS_TEST_CASE_IN_SUITE(GeometricalObjectsBinsSearchNearest, KratosCoreFastSuite) 
{
    constexpr double tolerance = 1e-12;

    Model current_model;

    // Cube coordinates
    const double cube_z = 0.3;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeSkinModelPart(current_model, 0.6, 0.9, cube_z);

    GeometricalObjectsBins bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd());

    double epsilon = 1.00e-6;
    Point near_point{epsilon,epsilon,epsilon};
    auto result = bins.SearchNearest(near_point);

    KRATOS_CHECK_NEAR(result.GetDistance(), (cube_z - epsilon), tolerance);

    std::size_t id = result.Get()->Id();
    KRATOS_CHECK(id == 3); 
}

/** Checks bins search nearest
*/
KRATOS_TEST_CASE_IN_SUITE(GeometricalObjectsBinsSearchNearestContainer, KratosCoreFastSuite) 
{
    constexpr double tolerance = 1e-12;

    Model current_model;
    
    // Cube coordinates
    const double cube_z = 0.3;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeSkinModelPart(current_model, 0.6, 0.9, cube_z);

    GeometricalObjectsBins bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd());

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    const double epsilon = 1.0e-6;
    auto p_node = r_point_model_part.CreateNewNode(1, epsilon, epsilon, epsilon);
    auto& r_array_nodes = r_point_model_part.Nodes();

    GeometricalObjectsBins::ResultTypeContainerMap results;
    bins.SearchNearest(r_array_nodes.begin(), r_array_nodes.end(), results);

    KRATOS_CHECK_EQUAL(results.NumberOfSearchResults(), 1);
    KRATOS_CHECK(results[*p_node].IsObjectFound());
    KRATOS_CHECK_EQUAL(results[*p_node].NumberOfGlobalResults(), 1);

    // Distances
    KRATOS_CHECK_NEAR(results[*p_node][0].GetDistance(), (cube_z - epsilon), tolerance);

    // Compute indices
    auto indices = results[*p_node].GetResultIndices();
    const std::size_t id = indices[0];
    KRATOS_CHECK_EQUAL(id, 3);
}

/** Checks bins empty search nearest 
*/
KRATOS_TEST_CASE_IN_SUITE(GeometricalObjectsBinsEmptySearchNearest, KratosCoreFastSuite) 
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = current_model.CreateModelPart("Skin");

    GeometricalObjectsBins bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd());

    Point center_point{0.00,0.00,0.00};
    auto result = bins.SearchNearest(center_point);

    KRATOS_CHECK_IS_FALSE(result.IsObjectFound());
}

/** Checks bins empty search nearest 
*/
KRATOS_TEST_CASE_IN_SUITE(GeometricalObjectsBinsEmptySearchNearestContainer, KratosCoreFastSuite) 
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = current_model.CreateModelPart("Skin");

    GeometricalObjectsBins bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd());

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    auto p_point = r_point_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto& r_array_nodes = r_point_model_part.Nodes();

    GeometricalObjectsBins::ResultTypeContainerMap results;
    bins.SearchNearest(r_array_nodes.begin(), r_array_nodes.end(), results);

    KRATOS_CHECK_EQUAL(results.NumberOfSearchResults(), 1);
    KRATOS_CHECK_IS_FALSE(results[*p_point].IsObjectFound());
}

/** Checks bins search is inside 
*/
KRATOS_TEST_CASE_IN_SUITE(GeometricalObjectsBinsSearchIsInside, KratosCoreFastSuite) 
{
    constexpr double tolerance = 1e-12;

    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = current_model.CreateModelPart("Skin");
    r_skin_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_skin_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_skin_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    r_skin_part.CreateNewNode(4, 0.0, 0.0, 1.0);
    Properties::Pointer p_properties(new Properties(0));
    r_skin_part.CreateNewElement("Element3D4N",  1, { 1,2,3,4 }, p_properties);

    GeometricalObjectsBins bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd());

    Point near_point{0.1,0.1,0.1};
    auto result = bins.SearchIsInside(near_point);

    KRATOS_CHECK(result.IsObjectFound());
    KRATOS_CHECK_NEAR(result.GetDistance(), 0.0, tolerance);
}

/** Checks bins search is inside 
*/
KRATOS_TEST_CASE_IN_SUITE(GeometricalObjectsBinsSearchIsInsideContainer, KratosCoreFastSuite) 
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeModelPart(current_model);

    GeometricalObjectsBins bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd());

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    auto p_inside_point = r_point_model_part.CreateNewNode(1, 0.5,0.5,0.5);
    auto& r_array_nodes = r_point_model_part.Nodes();

    GeometricalObjectsBins::ResultTypeContainerMap results;
    bins.SearchIsInside(r_array_nodes.begin(), r_array_nodes.end(), results);

    KRATOS_CHECK_EQUAL(results.NumberOfSearchResults(), 1);
    KRATOS_CHECK(results[*p_inside_point].IsObjectFound());
    KRATOS_CHECK_EQUAL(results[*p_inside_point].NumberOfGlobalResults(), 1);
}

/** Checks bins search is inside = not found
*/
KRATOS_TEST_CASE_IN_SUITE(GeometricalObjectsBinsSearchIsNotInside, KratosCoreFastSuite) 
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = current_model.CreateModelPart("Skin");
    r_skin_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_skin_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_skin_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    r_skin_part.CreateNewNode(4, 0.0, 0.0, 1.0);
    Properties::Pointer p_properties(new Properties(0));
    r_skin_part.CreateNewElement("Element3D4N",  1, { 1,2,3,4 }, p_properties);

    GeometricalObjectsBins bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd());

    Point near_point{0.5,0.5,0.5};
    auto result = bins.SearchIsInside(near_point);

    KRATOS_CHECK_IS_FALSE(result.IsObjectFound());
}

/** Checks bins search is inside = not found
*/
KRATOS_TEST_CASE_IN_SUITE(GeometricalObjectsBinsSearchIsNotInsideContainer, KratosCoreFastSuite) 
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeModelPart(current_model);

    GeometricalObjectsBins bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd());

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    auto p_outside_point = r_point_model_part.CreateNewNode(1, 100.0,100.0,100.0);
    auto& r_array_nodes = r_point_model_part.Nodes();

    GeometricalObjectsBins::ResultTypeContainerMap results;
    bins.SearchIsInside(r_array_nodes.begin(), r_array_nodes.end(), results);

    KRATOS_CHECK_EQUAL(results.NumberOfSearchResults(), 1);
    KRATOS_CHECK_IS_FALSE(results[*p_outside_point].IsObjectFound());
}

} // namespace Kratos::Testing.
