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

    KRATOS_CHECK_IS_FALSE(SearchUtilities::PointIsInsideBoundingBox(bounding_box, p_out_x));
    KRATOS_CHECK_IS_FALSE(SearchUtilities::PointIsInsideBoundingBox(bounding_box, p_out_y));
    KRATOS_CHECK_IS_FALSE(SearchUtilities::PointIsInsideBoundingBox(bounding_box, p_out_z));

    KRATOS_CHECK(SearchUtilities::PointIsInsideBoundingBox(bounding_box, p_in));
}

double GetBBoxValue(const int Index, const double Factor, const double Offset)
{
    return static_cast<double>(Index)*Factor - Offset;
}

KRATOS_TEST_CASE_IN_SUITE(SearchUtilitiesComputeBoundingBoxesWithTolerance, KratosCoreFastSuite)
{
    std::vector<double> bboxes_wrong_size(5);
    std::vector<double> bboxes_with_tol;

    KRATOS_DEBUG_CHECK_EXCEPTION_IS_THROWN(SearchUtilities::ComputeBoundingBoxesWithTolerance(bboxes_wrong_size, 1.235, bboxes_with_tol),
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
        KRATOS_CHECK_NEAR(bboxes_with_tol[i], (GetBBoxValue(i, factor, offset) + tolerance), 1e-12);

    for (int i=1; i<num_entries; i+=2)
        KRATOS_CHECK_NEAR(bboxes_with_tol[i], (GetBBoxValue(i, factor, offset) - tolerance), 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(SearchUtilitiesComputeBoundingBoxesWithToleranceCheckingNullBB, KratosCoreFastSuite)
{
    std::vector<double> bboxes_wrong_size(5);
    std::vector<double> bboxes_with_tol;

    KRATOS_DEBUG_CHECK_EXCEPTION_IS_THROWN(SearchUtilities::ComputeBoundingBoxesWithToleranceCheckingNullBB(bboxes_wrong_size, 1.235, bboxes_with_tol),
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
        KRATOS_CHECK_DOUBLE_EQUAL(bboxes_with_tol[i], 0.0);
    }
}

}  // namespace Kratos::Testing