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
#include "utilities/cpp_tests_utilities.h"

namespace Kratos::Testing 
{

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
    ModelPart& r_skin_part = CppTestsUtilities::CreateCubeSkinModelPart(current_model, cube_x, cube_y, cube_z);

    GeometricalObjectsBins bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd());

    auto bounding_box = bins.GetBoundingBox();

    KRATOS_EXPECT_NEAR(bounding_box.GetMinPoint()[0],-cube_x, tolerance);
    KRATOS_EXPECT_NEAR(bounding_box.GetMinPoint()[1],-cube_y, tolerance);
    KRATOS_EXPECT_NEAR(bounding_box.GetMinPoint()[2],-cube_z, tolerance);
    KRATOS_EXPECT_NEAR(bounding_box.GetMaxPoint()[0], cube_x, tolerance);
    KRATOS_EXPECT_NEAR(bounding_box.GetMaxPoint()[1], cube_y, tolerance);
    KRATOS_EXPECT_NEAR(bounding_box.GetMaxPoint()[2], cube_z, tolerance);
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
    ModelPart& r_skin_part = CppTestsUtilities::CreateCubeSkinModelPart(current_model, cube_x, cube_y, cube_z);

    GeometricalObjectsBins bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd());

    auto number_of_cells = bins.GetNumberOfCells();
    KRATOS_EXPECT_EQ(number_of_cells[0], 3);
    KRATOS_EXPECT_EQ(number_of_cells[1], 3);
    KRATOS_EXPECT_EQ(number_of_cells[2], 2);

    auto cell_sizes = bins.GetCellSizes();
    KRATOS_EXPECT_NEAR(cell_sizes[0], 2.00 * cube_x / 3.00, tolerance);
    KRATOS_EXPECT_NEAR(cell_sizes[1], 2.00 * cube_y / 3.00, tolerance);
    KRATOS_EXPECT_NEAR(cell_sizes[2], 2.00 * cube_z / 2.00, tolerance);
}


/** Checks bins AddObjectsToCells
*/
KRATOS_TEST_CASE_IN_SUITE(GeometricalObjectsBinsAddObjectsToCells, KratosCoreFastSuite) 
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = CppTestsUtilities::CreateCubeSkinModelPart(current_model);

    GeometricalObjectsBins bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd());

    auto& cell = bins.GetCell(0,0,0);
    KRATOS_EXPECT_EQ(cell.size(), 4);
    for(auto& geometrical_object : cell){
        std::size_t id = geometrical_object->GetId();
        KRATOS_EXPECT_TRUE((id == 1) ||(id == 2) ||(id == 7) ||(id == 11));
    }

    cell = bins.GetCell(2,2,1);
    KRATOS_EXPECT_EQ(cell.size(), 4);
    for(auto& geometrical_object : cell){
        std::size_t id = geometrical_object->GetId();
        KRATOS_EXPECT_TRUE((id == 3) ||(id == 4) ||(id == 6) ||(id == 10));
    }
}

/** Checks bins search in radius
*/
KRATOS_TEST_CASE_IN_SUITE(GeometricalObjectsBinsSearchInRadius, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = CppTestsUtilities::CreateCubeSkinModelPart(current_model);

    GeometricalObjectsBins bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd());

    std::vector<GeometricalObjectsBins::ResultType> results;
    Point center_point{0.00,0.00,0.00};

    bins.SearchInRadius(center_point, .29, results);
    KRATOS_EXPECT_EQ(results.size(), 0);

    bins.SearchInRadius(center_point, .3, results);
    KRATOS_EXPECT_EQ(results.size(), 4);

    bins.SearchInRadius(center_point, .4, results);
    KRATOS_EXPECT_EQ(results.size(), 4);

    bins.SearchInRadius(center_point, .6, results);
    KRATOS_EXPECT_EQ(results.size(), 8);

    bins.SearchInRadius(center_point, .7, results);
    KRATOS_EXPECT_EQ(results.size(), 8);

    bins.SearchInRadius(center_point, .9, results);
    KRATOS_EXPECT_EQ(results.size(), 12);
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
    ModelPart& r_skin_part = CppTestsUtilities::CreateCubeSkinModelPart(current_model, 0.6, 0.9, cube_z);

    GeometricalObjectsBins bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd());

    double epsilon = 1.00e-6;
    Point near_point{epsilon,epsilon,epsilon};
    auto result = bins.SearchNearestInRadius(near_point, cube_z - 1.e-4);

    KRATOS_EXPECT_FALSE(result.IsObjectFound());

    result = bins.SearchNearestInRadius(near_point, cube_z + 1.e-4);
    KRATOS_EXPECT_NEAR(result.GetDistance(), (cube_z - epsilon), tolerance);

    std::size_t id = result.Get()->Id();
    KRATOS_EXPECT_TRUE(id == 3);
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
    ModelPart& r_skin_part = CppTestsUtilities::CreateCubeSkinModelPart(current_model, 0.6, 0.9, cube_z);

    GeometricalObjectsBins bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd());

    double epsilon = 1.00e-6;
    Point near_point{epsilon,epsilon,epsilon};
    auto result = bins.SearchNearest(near_point);

    KRATOS_EXPECT_NEAR(result.GetDistance(), (cube_z - epsilon), tolerance);

    std::size_t id = result.Get()->Id();
    KRATOS_EXPECT_TRUE(id == 3); 
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

    KRATOS_EXPECT_FALSE(result.IsObjectFound());
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

    KRATOS_EXPECT_TRUE(result.IsObjectFound());
    KRATOS_EXPECT_NEAR(result.GetDistance(), 0.0, tolerance);
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

    KRATOS_EXPECT_FALSE(result.IsObjectFound());
}

} // namespace Kratos::Testing.
