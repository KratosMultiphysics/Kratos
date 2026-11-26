//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Juan I. Camarotti
//

#include "geometries/line_3d_2.h"
#include "testing/testing.h"
#include "custom_utilities/beam_mapper_utilities.h"
#include "custom_utilities/projection_utilities.h"
#include "mapping_application_variables.h"

namespace Kratos::Testing {

typedef Node NodeType;
typedef Geometry<NodeType> GeometryType;

namespace {

void TestHermitianShapeFunctionValues(
    const array_1d<double, 3>& local_coords,
    const std::vector<double>& exp_values,
    const std::vector<double>& exp_derivatives)
{
    Vector values, derivatives;
    BeamMapperUtilities::HermitianShapeFunctionsValues(values, derivatives, local_coords);

    KRATOS_EXPECT_EQ(values.size(), exp_values.size());
    KRATOS_EXPECT_EQ(derivatives.size(), exp_derivatives.size());

    for (std::size_t i = 0; i < exp_values.size(); ++i) {
        KRATOS_EXPECT_NEAR(exp_values[i], values[i], 1e-6);
        KRATOS_EXPECT_NEAR(exp_derivatives[i], derivatives[i], 1e-6);
    }
}

void TestProjectOnLineHermitian(
    GeometryType& rGeom,
    const Point& rPoint,
    const array_1d<double, 3>& local_coords_guess,
    const std::vector<double>& exp_values,
    const std::vector<double>& exp_derivatives,
    const double exp_proj_dist)
{
    Vector values, derivatives;
    Point proj_point;
    double proj_dist;
    ProjectionUtilities::PairingIndex pairing_index = 
        BeamMapperUtilities::ProjectOnLineHermitian(rGeom, rPoint, 1e-6, values, derivatives, proj_dist, proj_point);

    KRATOS_EXPECT_EQ(pairing_index, ProjectionUtilities::PairingIndex::Line_Inside);

    for (std::size_t i = 0; i < exp_values.size(); ++i) {
        KRATOS_EXPECT_NEAR(values[i], exp_values[i], 1e-6);
        KRATOS_EXPECT_NEAR(derivatives[i], exp_derivatives[i], 1e-6);
    }
    KRATOS_EXPECT_NEAR(proj_dist, exp_proj_dist, 1e-6);
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(BeamMapper_HermitianShapeFunctionsValues_Center, KratosMappingApplicationSerialTestSuite)
{
    array_1d<double, 3> local_coords;
    local_coords[0] = 0.5;
    local_coords[1] = 0.0;
    local_coords[2] = 0.0;

    // Correct Hermitian values and derivatives for ξ = 0.5
    std::vector<double> exp_values        {0.15625, 0.046875, 0.84375, -0.140625};
    std::vector<double> exp_derivatives   {-1.125, -0.3125, 1.125, 0.1875};

    TestHermitianShapeFunctionValues(local_coords, exp_values, exp_derivatives);
}

KRATOS_TEST_CASE_IN_SUITE(BeamMapper_ProjectOnLineHermitian_SimpleLine, KratosMappingApplicationSerialTestSuite)
{
    auto node_1(Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0));
    Line3D2<Node> line_geom(node_1, node_2);

    Point point_to_project(0.25, 0.0, 0.0);
    array_1d<double, 3> local_coords_guess = ZeroVector(3);

    std::vector<double> exp_values       {0.84375, 0.140625, 0.15625, -0.046875};
    std::vector<double> exp_derivatives  {-1.125, 0.1875, 1.125, -0.3125};
    double exp_proj_dist = 0.0;

    TestProjectOnLineHermitian(line_geom, point_to_project, local_coords_guess, exp_values, exp_derivatives, exp_proj_dist);
}

} // namespace Kratos::Testing
