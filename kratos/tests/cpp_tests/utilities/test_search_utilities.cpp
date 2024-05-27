//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//                   Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/point.h"
#include "utilities/search_utilities.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(SearchUtilitiesPointIsInsideBoundingBox, KratosCoreFastSuite)
{
    const SearchUtilities::BoundingBoxType bounding_box {10.5, -2.8, 3.89, -77.6, 4.64, 2.3};
    // xmax, xmin,  ymax, ymin,  zmax, zmin

    const Point p_out_x(10.6, 1.0, 3.8);
    const Point p_out_y(10.1, -80.0, 3.8);
    const Point p_out_z(10.1, 1.0, -3.8);
    const Point p_in(10.0, -30.78, 3.7);

    KRATOS_EXPECT_FALSE(SearchUtilities::PointIsInsideBoundingBox(bounding_box, p_out_x));
    KRATOS_EXPECT_FALSE(SearchUtilities::PointIsInsideBoundingBox(bounding_box, p_out_y));
    KRATOS_EXPECT_FALSE(SearchUtilities::PointIsInsideBoundingBox(bounding_box, p_out_z));

    KRATOS_EXPECT_TRUE(SearchUtilities::PointIsInsideBoundingBox(bounding_box, p_in));
}

double GetBBoxValue(const int Index, const double Factor, const double Offset)
{
    return static_cast<double>(Index)*Factor - Offset;
}

KRATOS_TEST_CASE_IN_SUITE(SearchUtilitiesComputeBoundingBoxesWithTolerance, KratosCoreFastSuite)
{
    std::vector<double> bboxes_wrong_size(5);
    std::vector<double> bboxes_with_tol;

    KRATOS_DEBUG_EXCEPT_EXCEPTION_IS_THROWN(SearchUtilities::ComputeBoundingBoxesWithTolerance(bboxes_wrong_size, 1.235, bboxes_with_tol),
        "Error: Bounding Boxes size has to be a multiple of 6!");

    // Cretae a vector containing the fake bboxes
    const int num_entries = 24;
    std::vector<double> bboxes(num_entries);

    const double factor = 1.2589;
    const double offset = 8.4;

    for (int i=0; i<num_entries; ++i)
        bboxes[i] = GetBBoxValue(i, factor, offset);

    const double tolerance = 5.478;

    SearchUtilities::ComputeBoundingBoxesWithTolerance(bboxes,
                                                       tolerance,
                                                       bboxes_with_tol);

    for (int i=0; i<num_entries; i+=2)
        KRATOS_EXPECT_NEAR(bboxes_with_tol[i], (GetBBoxValue(i, factor, offset) + tolerance), 1e-12);

    for (int i=1; i<num_entries; i+=2)
        KRATOS_EXPECT_NEAR(bboxes_with_tol[i], (GetBBoxValue(i, factor, offset) - tolerance), 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(SearchUtilitiesComputeBoundingBoxesWithToleranceCheckingNullBB, KratosCoreFastSuite)
{
    std::vector<double> bboxes_wrong_size(5);
    std::vector<double> bboxes_with_tol;

    KRATOS_DEBUG_EXCEPT_EXCEPTION_IS_THROWN(SearchUtilities::ComputeBoundingBoxesWithToleranceCheckingNullBB(bboxes_wrong_size, 1.235, bboxes_with_tol),
        "Error: Bounding Boxes size has to be a multiple of 6!");

    // Cretae a vector containing the fake bboxes
    const int num_entries = 24;
    std::vector<double> bboxes(num_entries, 0.0);

    const double tolerance = 5.478;

    SearchUtilities::ComputeBoundingBoxesWithToleranceCheckingNullBB(bboxes,
                                                       tolerance,
                                                       bboxes_with_tol);
    // Check that the bboxes are all zero
    for (int i=0; i<num_entries; ++i) {
        KRATOS_EXPECT_DOUBLE_EQ(bboxes_with_tol[i], 0.0);
    }
}

KRATOS_TEST_CASE_IN_SUITE(SynchronousPointSynchronization, KratosCoreFastSuite) 
{
    // The data communicator
    const DataCommunicator& r_data_comm = DataCommunicator();

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

KRATOS_TEST_CASE_IN_SUITE(SynchronousPointSynchronizationWithRadius, KratosCoreFastSuite) 
{
    // The data communicator
    const DataCommunicator& r_data_comm = DataCommunicator();

    // Data communicator
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

}  // namespace Kratos::Testing