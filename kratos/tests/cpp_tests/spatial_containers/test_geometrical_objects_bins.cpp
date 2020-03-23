//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
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

namespace Kratos {
namespace Testing {

    ModelPart& CreateCubeSkinModelPart(Model& CurrentModel, double HalfX, double HalfY, double HalfZ){
        // Generate the cube skin
        ModelPart& skin_part = CurrentModel.CreateModelPart("Skin");
        skin_part.CreateNewNode(1, -HalfX, -HalfY, -HalfZ);
        skin_part.CreateNewNode(2,  HalfX, -HalfY, -HalfZ);
        skin_part.CreateNewNode(3,  HalfX,  HalfY, -HalfZ);
        skin_part.CreateNewNode(4, -HalfX,  HalfY, -HalfZ);
        skin_part.CreateNewNode(5, -HalfX, -HalfY,  HalfZ);
        skin_part.CreateNewNode(6,  HalfX, -HalfY,  HalfZ);
        skin_part.CreateNewNode(7,  HalfX,  HalfY,  HalfZ);
        skin_part.CreateNewNode(8, -HalfX,  HalfY,  HalfZ);
        Properties::Pointer p_properties(new Properties(0));
        skin_part.CreateNewElement("Element3D3N",  1, { 1,2,3 }, p_properties);
        skin_part.CreateNewElement("Element3D3N",  2, { 1,3,4 }, p_properties);
        skin_part.CreateNewElement("Element3D3N",  3, { 5,6,7 }, p_properties);
        skin_part.CreateNewElement("Element3D3N",  4, { 5,7,8 }, p_properties);
        skin_part.CreateNewElement("Element3D3N",  5, { 3,6,2 }, p_properties);
        skin_part.CreateNewElement("Element3D3N",  6, { 3,7,6 }, p_properties);
        skin_part.CreateNewElement("Element3D3N",  7, { 4,5,1 }, p_properties);
        skin_part.CreateNewElement("Element3D3N",  8, { 4,8,5 }, p_properties);
        skin_part.CreateNewElement("Element3D3N",  9, { 3,4,8 }, p_properties);
        skin_part.CreateNewElement("Element3D3N", 10, { 3,8,7 }, p_properties);
        skin_part.CreateNewElement("Element3D3N", 11, { 2,1,5 }, p_properties);
        skin_part.CreateNewElement("Element3D3N", 12, { 2,5,6 }, p_properties);

        return skin_part;
    }

    /** Checks bins bounding box
    */
    KRATOS_TEST_CASE_IN_SUITE(GeometricalObjectsBinsBoundingBox, KratosFastSuite) {
        constexpr double tolerance = 1e-12;

        Model current_model;

        const double cube_x = 0.6;
        const double cube_y = 0.9;
        const double cube_z = 0.3;

        // Generate the cube skin
        ModelPart& skin_part = CreateCubeSkinModelPart(current_model, cube_x, cube_y, cube_z);

        GeometricalObjectsBins bins(skin_part.ElementsBegin(), skin_part.ElementsEnd());

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
    KRATOS_TEST_CASE_IN_SUITE(GeometricalObjectsBinsCellSizes, KratosFastSuite) {
        constexpr double tolerance = 1e-12;

        Model current_model;

        const double cube_x = 0.6;
        const double cube_y = 0.9;
        const double cube_z = 0.3;

        // Generate the cube skin
        ModelPart& skin_part = CreateCubeSkinModelPart(current_model, cube_x, cube_y, cube_z);

        GeometricalObjectsBins bins(skin_part.ElementsBegin(), skin_part.ElementsEnd());

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
    KRATOS_TEST_CASE_IN_SUITE(GeometricalObjectsBinsAddObjectsToCells, KratosFastSuite) {
        Model current_model;

        const double cube_x = 0.6;
        const double cube_y = 0.9;
        const double cube_z = 0.3;

        // Generate the cube skin
        ModelPart& skin_part = CreateCubeSkinModelPart(current_model, cube_x, cube_y, cube_z);

        GeometricalObjectsBins bins(skin_part.ElementsBegin(), skin_part.ElementsEnd());

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
    KRATOS_TEST_CASE_IN_SUITE(GeometricalObjectsBinsSearchInRadius, KratosFastSuite) {
        Model current_model;

        const double cube_x = 0.6;
        const double cube_y = 0.9;
        const double cube_z = 0.3;

        // Generate the cube skin
        ModelPart& skin_part = CreateCubeSkinModelPart(current_model, cube_x, cube_y, cube_z);

        GeometricalObjectsBins bins(skin_part.ElementsBegin(), skin_part.ElementsEnd());

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
 
    /** Checks bins search nearest
    */
    KRATOS_TEST_CASE_IN_SUITE(GeometricalObjectsBinsSearchNearestInRadius, KratosFastSuite) {
        constexpr double tolerance = 1e-12;

        Model current_model;

        const double cube_x = 0.6;
        const double cube_y = 0.9;
        const double cube_z = 0.3;

        // Generate the cube skin
        ModelPart& skin_part = CreateCubeSkinModelPart(current_model, cube_x, cube_y, cube_z);

        GeometricalObjectsBins bins(skin_part.ElementsBegin(), skin_part.ElementsEnd());

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
    KRATOS_TEST_CASE_IN_SUITE(GeometricalObjectsBinsSearchNearest, KratosFastSuite) {
        constexpr double tolerance = 1e-12;

        Model current_model;

        const double cube_x = 0.6;
        const double cube_y = 0.9;
        const double cube_z = 0.3;

        // Generate the cube skin
        ModelPart& skin_part = CreateCubeSkinModelPart(current_model, cube_x, cube_y, cube_z);

        GeometricalObjectsBins bins(skin_part.ElementsBegin(), skin_part.ElementsEnd());

        double epsilon = 1.00e-6;
        Point near_point{epsilon,epsilon,epsilon};
        auto result = bins.SearchNearest(near_point);

        KRATOS_CHECK_NEAR(result.GetDistance(), (cube_z - epsilon), tolerance);
        
        std::size_t id = result.Get()->Id();
        KRATOS_CHECK(id == 3); 
   }

    /** Checks bins empty search nearest 
    */
    KRATOS_TEST_CASE_IN_SUITE(GeometricalObjectsBinsEmptySearchNearest, KratosFastSuite) {
        Model current_model;

        // Generate the cube skin
        ModelPart& skin_part = current_model.CreateModelPart("Skin");

        GeometricalObjectsBins bins(skin_part.ElementsBegin(), skin_part.ElementsEnd());

        Point center_point{0.00,0.00,0.00};
        auto result = bins.SearchNearest(center_point);

        KRATOS_CHECK_IS_FALSE(result.IsObjectFound());
    }

} // namespace Testing.
} // namespace Kratos.
