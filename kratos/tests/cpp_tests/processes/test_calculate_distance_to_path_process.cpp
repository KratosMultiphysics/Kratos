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

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/checks.h"
#include "geometries/line_3d_2.h"
#include "utilities/geometrical_projection_utilities.h"
#include "processes/calculate_distance_to_path_process.h"

namespace Kratos::Testing
{
    constexpr double TOLERANCE_DISTANCE_PATH = 1e-6;

    KRATOS_TEST_CASE_IN_SUITE(FastMinimalDistanceOnLineWithRadius, KratosCoreFastSuite)
	{
        double distance;
        double radius = 0.0;
        auto line = Kratos::make_shared<Line3D2<Node<3>>>(
            Kratos::make_intrusive<Node<3>>(1, 0.0, 0.0, 0.0),
            Kratos::make_intrusive<Node<3>>(2, 1.0, 0.0, 0.0)
        );
        Point point1(0.0,0.0,0.1);
        distance = GeometricalProjectionUtilities::FastMinimalDistanceOnLine(*line, point1);
        auto distance_computed_type = CalculateDistanceToPathProcess<true>::FastMinimalDistanceOnLineWithRadius(distance, *line, point1, radius);
        KRATOS_CHECK_NEAR(distance, 0.1, TOLERANCE_DISTANCE_PATH);
        KRATOS_CHECK_EQUAL(distance_computed_type, DistanceComputed::NO_RADIUS);

        radius = 0.01;
        Point point2(0.0,0.0,0.1);
        distance = GeometricalProjectionUtilities::FastMinimalDistanceOnLine(*line, point2);
        distance_computed_type = CalculateDistanceToPathProcess<true>::FastMinimalDistanceOnLineWithRadius(distance, *line, point2, radius);
        KRATOS_CHECK_NEAR(distance, 0.09, TOLERANCE_DISTANCE_PATH);
        KRATOS_CHECK_EQUAL(distance_computed_type, DistanceComputed::RADIUS_PROJECTED);

        Point point3(-0.1,0.0,0.1);
        distance = GeometricalProjectionUtilities::FastMinimalDistanceOnLine(*line, point3);
        distance_computed_type = CalculateDistanceToPathProcess<true>::FastMinimalDistanceOnLineWithRadius(distance, *line, point3, radius);
        KRATOS_CHECK_NEAR(distance, std::sqrt(std::pow(0.09, 2) * 2), TOLERANCE_DISTANCE_PATH);
        KRATOS_CHECK_EQUAL(distance_computed_type, DistanceComputed::RADIUS_NOT_PROJECTED_OUTSIDE);

        radius = 0.1;
        Point point4(-0.1,0.0,0.09);
        distance = GeometricalProjectionUtilities::FastMinimalDistanceOnLine(*line, point4);
        distance_computed_type = CalculateDistanceToPathProcess<true>::FastMinimalDistanceOnLineWithRadius(distance, *line, point4, radius);
        KRATOS_CHECK_NEAR(distance, -(std::sqrt(std::pow(0.01, 2) + std::pow(0.1, 2))), TOLERANCE_DISTANCE_PATH);
        KRATOS_CHECK_EQUAL(distance_computed_type, DistanceComputed::RADIUS_NOT_PROJECTED_INSIDE);
    }
}