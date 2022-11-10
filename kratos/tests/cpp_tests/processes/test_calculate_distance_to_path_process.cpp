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
    constexpr double TOLERANCE = 1e-6;

    KRATOS_TEST_CASE_IN_SUITE(FastMinimalDistanceOnLineWithRadius, KratosCoreFastSuite)
	{
        double distance;
        double radius = 0.0;
        auto line = Kratos::make_shared<Line3D2<Node<3>>>(
            Kratos::make_intrusive<Node<3>>(0.0, 0.0, 0.0),
            Kratos::make_intrusive<Node<3>>(1.0, 0.0, 0.0)
        );
        Point point(0.0,0.0,0.1);
        distance = GeometricalProjectionUtilities::FastMinimalDistanceOnLine(*line, point);
        const auto distance_computed_type = CalculateDistanceToPathProcess<true>::FastMinimalDistanceOnLineWithRadius(distance, *line, point, radius);
        KRATOS_CHECK_NEAR(distance, 0.1, TOLERANCE);
    }
}