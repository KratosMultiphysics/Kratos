// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
// System includes

// External includes

// Project includes
#include "containers/array_1d.h"

// Application includes
#include "custom_utilities/metrics_math_utils.h"
#include "tests/cpp_tests/meshing_fast_suite.h"

namespace Kratos::Testing 
{

static constexpr double TOLERANCE = 1.0e-5;

/**
 * Checks the correct work of the IntersectMetrics
 * Test 2D
 */
KRATOS_TEST_CASE_IN_SUITE(IntersectMetrics2D, KratosMeshingApplicationFastSuite)
{
    array_1d<double, 3> metric_1, metric_2, intersected_metric;

    metric_1[0] = 1.0;
    metric_1[1] = 2.0;
    metric_1[2] = 0.0;

    metric_2[0] = 2.0;
    metric_2[1] = 1.0;
    metric_2[2] = 0.0;

    noalias(intersected_metric) = MetricsMathUtils<2>::IntersectMetrics(metric_1, metric_2);

    KRATOS_EXPECT_NEAR(intersected_metric[0], 2.0, TOLERANCE);
    KRATOS_EXPECT_NEAR(intersected_metric[1], 2.0, TOLERANCE);
    KRATOS_EXPECT_NEAR(intersected_metric[2], 0.0, TOLERANCE);

    metric_1[2] = 0.1;

    metric_2[2] = 0.2;

    noalias(intersected_metric) = MetricsMathUtils<2>::IntersectMetrics(metric_1, metric_2);

    KRATOS_EXPECT_NEAR(intersected_metric[0], 1.95164, TOLERANCE);
    KRATOS_EXPECT_NEAR(intersected_metric[1], 2.00933, TOLERANCE);
    KRATOS_EXPECT_NEAR(intersected_metric[2], 0.00576923, TOLERANCE);
}

/**
 * Checks the correct work of the IntersectMetrics
 * Test 3D
 */
KRATOS_TEST_CASE_IN_SUITE(IntersectMetrics3D, KratosMeshingApplicationFastSuite)
{
    array_1d<double, 6> metric_1, metric_2, intersected_metric;

    metric_1[0] = 1.0;
    metric_1[1] = 2.0;
    metric_1[2] = 2.0;
    metric_1[3] = 0.0;
    metric_1[4] = 0.0;
    metric_1[5] = 0.0;

    metric_2[0] = 2.0;
    metric_2[1] = 1.0;
    metric_2[2] = 1.0;
    metric_2[3] = 0.0;
    metric_2[4] = 0.0;
    metric_2[5] = 0.0;

    noalias(intersected_metric) = MetricsMathUtils<3>::IntersectMetrics(metric_1, metric_2);

    KRATOS_EXPECT_NEAR(intersected_metric[0], 2.0, TOLERANCE);
    KRATOS_EXPECT_NEAR(intersected_metric[1], 2.0, TOLERANCE);
    KRATOS_EXPECT_NEAR(intersected_metric[2], 2.0, TOLERANCE);
    KRATOS_EXPECT_NEAR(intersected_metric[3], 0.0, TOLERANCE);
    KRATOS_EXPECT_NEAR(intersected_metric[4], 0.0, TOLERANCE);
    KRATOS_EXPECT_NEAR(intersected_metric[5], 0.0, TOLERANCE);

    metric_1[3] = 0.1;
    metric_1[4] = 0.3;
    metric_1[5] = 0.1;

    metric_2[3] = 0.2;
    metric_2[4] = 0.1;
    metric_2[5] = 0.2;

    noalias(intersected_metric) = MetricsMathUtils<3>::IntersectMetrics(metric_1, metric_2);

    KRATOS_EXPECT_NEAR(intersected_metric[0], 1.92815, TOLERANCE);
    KRATOS_EXPECT_NEAR(intersected_metric[1], 2.00271, TOLERANCE);
    KRATOS_EXPECT_NEAR(intersected_metric[2], 2.00417, TOLERANCE);
    KRATOS_EXPECT_NEAR(intersected_metric[3], 0.0105435, TOLERANCE);
    KRATOS_EXPECT_NEAR(intersected_metric[4], 0.307946, TOLERANCE);
    KRATOS_EXPECT_NEAR(intersected_metric[5], 0.0427123, TOLERANCE);
}

} // namespace Kratos::Testing
