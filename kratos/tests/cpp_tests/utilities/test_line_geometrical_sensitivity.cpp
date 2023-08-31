//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:
//

#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"
#include "geometries/triangle_3d_3.h"
#include "includes/expect.h"
#include "utilities/line_sensitivity_utility.h"
#include "utilities/math_utils.h"
#include "testing/testing.h"

namespace Kratos
{
namespace Testing
{

Geometry<Point>::Pointer CreateLine2D2N()
{
    Geometry<Point>::PointsArrayType points;
    points.push_back(Kratos::make_shared<Point>(-1.2, 0.5, 0.0));
    points.push_back(Kratos::make_shared<Point>(3.1, -1.0, 0.0));

    return Geometry<Point>::Pointer(new Line2D2<Point>(points));
}

Geometry<Point>::Pointer CreateLine2D3N()
{
    Geometry<Point>::PointsArrayType points;
    points.push_back(Kratos::make_shared<Point>(-1.2,  0.5, 0.0));
    points.push_back(Kratos::make_shared<Point>( 0.7,  0.6, 0.0));
    points.push_back(Kratos::make_shared<Point>( 3.1, -1.0, 0.0));

    return Geometry<Point>::Pointer(new Line2D3<Point>(points));
}

Geometry<Point>::Pointer CreateLine3D2N()
{
    Geometry<Point>::PointsArrayType points;
    points.push_back(Kratos::make_shared<Point>(-1.2, 0.5, -1.1));
    points.push_back(Kratos::make_shared<Point>(3.1, -1.0, 0.9));

    return Geometry<Point>::Pointer(new Line3D2<Point>(points));
}

Geometry<Point>::Pointer CreateLine3D3N()
{
    Geometry<Point>::PointsArrayType points;
    points.push_back(Kratos::make_shared<Point>(-1.2,  0.5, -1.1));
    points.push_back(Kratos::make_shared<Point>( 0.7,  0.6,  1.3));
    points.push_back(Kratos::make_shared<Point>( 3.1, -1.0,  0.9));

    return Geometry<Point>::Pointer(new Line3D3<Point>(points));
}

Geometry<Point>::Pointer CreateTriangle3D3N()
{
    Geometry<Point>::PointsArrayType points;
    points.push_back(Kratos::make_shared<Point>(-1.2,  0.5, -1.1));
    points.push_back(Kratos::make_shared<Point>( 0.7,  0.6,  1.3));
    points.push_back(Kratos::make_shared<Point>( 3.1, -1.0,  0.9));

    return Geometry<Point>::Pointer(new Triangle3D3<Point>(points));
}

DenseMatrix<double> GetJacobian(
    const Geometry<Point>& rGeometry,
    GeometryData::IntegrationMethod Quadrature,
    unsigned int GaussPointIndex)
{
    const auto& rDN_De = rGeometry.ShapeFunctionLocalGradient(GaussPointIndex, Quadrature);
    DenseMatrix<double> jacobian(rGeometry.WorkingSpaceDimension(), rGeometry.LocalSpaceDimension());
    DenseMatrix<double> coordinates(rGeometry.WorkingSpaceDimension(), rGeometry.PointsNumber());

    for (unsigned int i = 0; i < rGeometry.PointsNumber(); i++)
    {
        const auto& r_coordinates = rGeometry[i].Coordinates();
        for (unsigned int d = 0; d < rGeometry.WorkingSpaceDimension(); d++)
        {
            coordinates(d,i) = r_coordinates[d];
        }
    }

    noalias(jacobian) = prod(coordinates,rDN_De);
    return jacobian;
}

double IntegrationPointWeight(
    const Geometry<Point>& rGeometry,
    GeometryData::IntegrationMethod Quadrature,
    unsigned int GaussPointIndex)
{
    auto jacobian = GetJacobian(rGeometry, Quadrature, GaussPointIndex);
    DenseMatrix<double> aux = prod(trans(jacobian), jacobian);
    double determinant = MathUtils<double>::Det(aux);
    return rGeometry.IntegrationPoints(Quadrature)[GaussPointIndex].Weight() * std::sqrt(determinant);
}

void CheckIntegrationPointWeightSensitivity(
    Geometry<Point>::Pointer pGeometry,
    GeometryData::IntegrationMethod Quadrature,
    double Perturbation = 1e-7,
    double Tolerance = 1e-7)
{
    auto& r_geometry = *pGeometry;
    const auto& integration_points = r_geometry.IntegrationPoints(Quadrature);
    const unsigned int num_integration_points = integration_points.size();

    for (unsigned int g = 0; g < num_integration_points; g++)
    {
        const auto& r_DN_De = r_geometry.ShapeFunctionLocalGradient(g, Quadrature);
        const auto jacobian = GetJacobian(r_geometry, Quadrature, g);
        LineSensitivityUtility line_sensitivity(jacobian, r_DN_De);
        double result;

        double base_weight = IntegrationPointWeight(r_geometry, Quadrature, g);
        ShapeParameter::Sequence s(r_geometry.PointsNumber(), r_geometry.WorkingSpaceDimension());
        while(s)
        {
            const auto& deriv = s.CurrentValue();
            r_geometry[deriv.NodeIndex].Coordinates()[deriv.Direction] += Perturbation;
            double perturbed_weight = IntegrationPointWeight(r_geometry, Quadrature, g);

            double finite_difference_sensitivity = (perturbed_weight-base_weight)/Perturbation;
            line_sensitivity.CalculateSensitivity(deriv, result);
            double obtained_sensitivity = result * integration_points[g].Weight();

            KRATOS_EXPECT_NEAR(obtained_sensitivity, finite_difference_sensitivity, Tolerance);

            // undo perturbation for next step
            r_geometry[deriv.NodeIndex].Coordinates()[deriv.Direction] -= Perturbation;
            ++s;
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(LineGeometricalSensitivity_Line2D2N, KratosCoreFastSuite)
{
    auto p_geometry = CreateLine2D2N();
    CheckIntegrationPointWeightSensitivity(p_geometry, GeometryData::IntegrationMethod::GI_GAUSS_1);
    CheckIntegrationPointWeightSensitivity(p_geometry, GeometryData::IntegrationMethod::GI_GAUSS_2);
    CheckIntegrationPointWeightSensitivity(p_geometry, GeometryData::IntegrationMethod::GI_GAUSS_3);
    CheckIntegrationPointWeightSensitivity(p_geometry, GeometryData::IntegrationMethod::GI_GAUSS_4);
}

KRATOS_TEST_CASE_IN_SUITE(LineGeometricalSensitivity_Line2D3N, KratosCoreFastSuite)
{
    auto p_geometry = CreateLine2D3N();
    CheckIntegrationPointWeightSensitivity(p_geometry, GeometryData::IntegrationMethod::GI_GAUSS_1);
    CheckIntegrationPointWeightSensitivity(p_geometry, GeometryData::IntegrationMethod::GI_GAUSS_2);
    CheckIntegrationPointWeightSensitivity(p_geometry, GeometryData::IntegrationMethod::GI_GAUSS_3);
    CheckIntegrationPointWeightSensitivity(p_geometry, GeometryData::IntegrationMethod::GI_GAUSS_4);
}

KRATOS_TEST_CASE_IN_SUITE(LineGeometricalSensitivity_Line3D2N, KratosCoreFastSuite)
{
    auto p_geometry = CreateLine3D2N();
    CheckIntegrationPointWeightSensitivity(p_geometry, GeometryData::IntegrationMethod::GI_GAUSS_1);
    CheckIntegrationPointWeightSensitivity(p_geometry, GeometryData::IntegrationMethod::GI_GAUSS_2);
    CheckIntegrationPointWeightSensitivity(p_geometry, GeometryData::IntegrationMethod::GI_GAUSS_3);
    CheckIntegrationPointWeightSensitivity(p_geometry, GeometryData::IntegrationMethod::GI_GAUSS_4);
}

KRATOS_TEST_CASE_IN_SUITE(LineGeometricalSensitivity_Line3D3N, KratosCoreFastSuite)
{
    auto p_geometry = CreateLine3D3N();
    CheckIntegrationPointWeightSensitivity(p_geometry, GeometryData::IntegrationMethod::GI_GAUSS_1);
    CheckIntegrationPointWeightSensitivity(p_geometry, GeometryData::IntegrationMethod::GI_GAUSS_2);
    CheckIntegrationPointWeightSensitivity(p_geometry, GeometryData::IntegrationMethod::GI_GAUSS_3);
    CheckIntegrationPointWeightSensitivity(p_geometry, GeometryData::IntegrationMethod::GI_GAUSS_4);
}

KRATOS_TEST_CASE_IN_SUITE(SurfaceGeometricalSensitivity_Triangle3D3N, KratosCoreFastSuite)
{
    auto p_geometry = CreateTriangle3D3N();
    CheckIntegrationPointWeightSensitivity(p_geometry, GeometryData::IntegrationMethod::GI_GAUSS_1);
    CheckIntegrationPointWeightSensitivity(p_geometry, GeometryData::IntegrationMethod::GI_GAUSS_2);
    CheckIntegrationPointWeightSensitivity(p_geometry, GeometryData::IntegrationMethod::GI_GAUSS_3);
    CheckIntegrationPointWeightSensitivity(p_geometry, GeometryData::IntegrationMethod::GI_GAUSS_4);
}

} // namespace Testing
} // namespace Kratos
