//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/line_2d_2.h"
#include "spatial_containers/spatial_search_result.h"
#include "spatial_containers/spatial_search_result_container_map.h"

namespace Kratos::Testing 
{

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerMapInitializeResult, KratosMPICoreFastSuite)
{    
    // Create a test object
    SpatialSearchResultContainerMap<GeometricalObject> container_map;

    // Initialize result
    Point point = Point(0.5, 0.0, 0.0);
    container_map.InitializeResult(point);

    // Check that the result was added correctly
    KRATOS_EXPECT_TRUE(container_map.HasResult(point));
    Point fake_point = Point(1.5, 0.0, 0.0);
    KRATOS_EXPECT_FALSE(container_map.HasResult(fake_point));
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerMapClear, KratosMPICoreFastSuite)
{    
    // Create a test object
    SpatialSearchResultContainerMap<GeometricalObject> container_map;

    // Initialize result
    Point point = Point(0.5, 0.0, 0.0);
    container_map.InitializeResult(point);

    // Check that the result was added correctly
    KRATOS_EXPECT_TRUE(container_map.HasResult(point));
    container_map.Clear();
    KRATOS_EXPECT_FALSE(container_map.HasResult(point));
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerMapOperators, KratosMPICoreFastSuite)
{    
    // Create a test object
    SpatialSearchResultContainerMap<GeometricalObject> container_map;

    // Initialize result
    Point point = Point(0.5, 0.0, 0.0);
    container_map.InitializeResult(point);

    // Check that the result was added correctly
    auto& r_result = container_map[point];
    auto& r_local_pointers = r_result.GetLocalResults();
    KRATOS_EXPECT_EQ(r_local_pointers.size(), 0);
    KRATOS_EXPECT_EQ(r_local_pointers.size(), r_result.NumberOfLocalResults());
}

}  // namespace Kratos::Testing