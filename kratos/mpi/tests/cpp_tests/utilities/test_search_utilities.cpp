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
#include "geometries/point.h"
#include "utilities/search_utilities.h"

namespace Kratos::Testing 
{

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISynchronousPointSynchronization, KratosMPICoreFastSuite) 
{
    // The data communicator
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    // MPI data communicator
    const int rank = r_data_comm.Rank();
    const int world_size = r_data_comm.Size();

    // Define test data
    double value = static_cast<double>(rank);
    double value2 = 2 * value;
    std::vector<Point> points({Point(value, value, value), Point(value2, value2, value2)});

    // Call the function
    std::vector<double> all_points_coordinates;
    std::vector<long int> all_points_ids;
    SearchUtilities::SynchronousPointSynchronization(points.begin(), points.end(), all_points_coordinates, all_points_ids, r_data_comm);

    // Check the results
    KRATOS_EXPECT_EQ(static_cast<int>(all_points_coordinates.size()), 2 * 3 * world_size);
    KRATOS_EXPECT_EQ(static_cast<int>(all_points_ids.size()), 2 * world_size);
    for (int i_rank = 0; i_rank < world_size; ++i_rank) {
        value = static_cast<double>(i_rank);
        value2 = 2 * value;
        KRATOS_EXPECT_EQ(all_points_ids[i_rank * 2    ], i_rank * 2);
        KRATOS_EXPECT_EQ(all_points_ids[i_rank * 2 + 1], i_rank * 2 + 1);
        for (int j = 0; j < 3; ++j) {
            KRATOS_EXPECT_DOUBLE_EQ(all_points_coordinates[i_rank * 6 + j    ], value );
            KRATOS_EXPECT_DOUBLE_EQ(all_points_coordinates[i_rank * 6 + j + 3], value2);
        }
    }
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISynchronousPointSynchronizationWithRadius, KratosMPICoreFastSuite) 
{
    // The data communicator
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    // MPI data communicator
    const int rank = r_data_comm.Rank();
    const int world_size = r_data_comm.Size();

    // Define test data
    double value = static_cast<double>(rank);
    double value2 = 2 * value;
    std::vector<Point> points({Point(value, value, value), Point(value2, value2, value2)});
    std::vector<double> local_radius = {value, value2};

    // Call the function
    std::vector<double> all_points_coordinates;
    std::vector<long int> all_points_ids;
    auto radius = SearchUtilities::SynchronousPointSynchronizationWithRadius(points.begin(), points.end(), all_points_coordinates, all_points_ids, local_radius, r_data_comm);

    // Check the results
    KRATOS_EXPECT_EQ(static_cast<int>(all_points_coordinates.size()), 2 * 3 * world_size);
    KRATOS_EXPECT_EQ(static_cast<int>(all_points_ids.size()), 2 * world_size);
    for (int i_rank = 0; i_rank < world_size; ++i_rank) {
        value = static_cast<double>(i_rank);
        value2 = 2 * value;
        KRATOS_EXPECT_DOUBLE_EQ(radius[i_rank * 2    ], value );
        KRATOS_EXPECT_DOUBLE_EQ(radius[i_rank * 2 + 1], value2);
        KRATOS_EXPECT_EQ(all_points_ids[i_rank * 2    ], i_rank * 2);
        KRATOS_EXPECT_EQ(all_points_ids[i_rank * 2 + 1], i_rank * 2 + 1);
        for (int j = 0; j < 3; ++j) {
            KRATOS_EXPECT_DOUBLE_EQ(all_points_coordinates[i_rank * 6 + j    ], value );
            KRATOS_EXPECT_DOUBLE_EQ(all_points_coordinates[i_rank * 6 + j + 3], value2);
        }
    }
}

}  // namespace Kratos::Testing.