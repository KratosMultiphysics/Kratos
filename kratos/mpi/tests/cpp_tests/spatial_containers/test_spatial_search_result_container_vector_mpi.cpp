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
#include "spatial_containers/spatial_search_result_container_vector.h"

namespace Kratos::Testing 
{

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerVectorInitializeResult, KratosMPICoreFastSuite)
{    
    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_map;

    // Initialize result
    const std::size_t index = 0;
    container_map.InitializeResult(index);

    // Check that the result was added correctly
    KRATOS_EXPECT_TRUE(container_map.HasResult(index));
    const std::size_t fake_index = 1;
    KRATOS_EXPECT_FALSE(container_map.HasResult(fake_index));
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerVectorClear, KratosMPICoreFastSuite)
{    
    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_map;

    // Initialize result
    Point point = Point(0.5, 0.0, 0.0);
    const std::size_t index = 0;
    container_map.InitializeResult(index);

    // Check that the result was added correctly
    KRATOS_EXPECT_TRUE(container_map.HasResult(index));
    container_map.Clear();
    KRATOS_EXPECT_FALSE(container_map.HasResult(index));
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerVectorOperators, KratosMPICoreFastSuite)
{    
    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_map;

    // Initialize result
    const std::size_t index = 0;
    container_map.InitializeResult(index);

    // Check that the result was added correctly
    auto& r_result = container_map[index];
    auto& r_local_pointers = r_result.GetLocalResults();
    KRATOS_EXPECT_EQ(r_local_pointers.size(), 0);
    KRATOS_EXPECT_EQ(r_local_pointers.size(), r_result.NumberOfLocalResults());
}

}  // namespace Kratos::Testing