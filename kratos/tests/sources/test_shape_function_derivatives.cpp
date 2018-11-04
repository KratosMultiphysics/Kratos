#include "tests/geometries/test_shape_function_derivatives.h"

#include <functional>
#include "containers/array_1d.h"
#include "testing/testing.h"
#include "tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

void TestShapeFunctionsLocalGradients(Geometry<Node<3>> const& rGeom)
{
    KRATOS_TRY;
    auto local_gradients = rGeom.ShapeFunctionsLocalGradients();
    auto integration_points = rGeom.IntegrationPoints();
    for (unsigned g = 0; g < rGeom.IntegrationPointsNumber(); ++g)
        TestShapeFunctionsLocalGradient(rGeom, integration_points[g], local_gradients[g]);
    KRATOS_CATCH("");
}

void TestShapeFunctionsLocalGradients(Geometry<Node<3>> const& rGeom,
                                      GeometryData::IntegrationMethod ThisMethod)
{
    KRATOS_TRY;
    auto local_gradients = rGeom.ShapeFunctionsLocalGradients(ThisMethod);
    auto integration_points = rGeom.IntegrationPoints(ThisMethod);
    for (unsigned g = 0; g < rGeom.IntegrationPointsNumber(ThisMethod); ++g)
        TestShapeFunctionsLocalGradient(rGeom, integration_points[g], local_gradients[g]);
    KRATOS_CATCH("");
}

void TestShapeFunctionsLocalGradients_IntegrationPointIndex(Geometry<Node<3>> const& rGeom)
{
    KRATOS_TRY;
    auto integration_points = rGeom.IntegrationPoints();
    for (std::size_t i = 0; i < integration_points.size(); ++i)
        TestShapeFunctionsLocalGradient(rGeom, integration_points[i],
                                        rGeom.ShapeFunctionLocalGradient(i));
    KRATOS_CATCH("");
}

void TestShapeFunctionsLocalGradients_IntegrationPointIndex(
    Geometry<Node<3>> const& rGeom, GeometryData::IntegrationMethod ThisMethod)
{
    KRATOS_TRY;
    auto integration_points = rGeom.IntegrationPoints(ThisMethod);
    for (std::size_t i = 0; i < integration_points.size(); ++i)
        TestShapeFunctionsLocalGradient(rGeom, integration_points[i],
                                        rGeom.ShapeFunctionLocalGradient(i, ThisMethod));
    KRATOS_CATCH("");
}


void TestShapeFunctionsLocalGradients_IntegrationPointIndex_ShapeFunctionIndex(
    Geometry<Node<3>> const& rGeom, GeometryData::IntegrationMethod ThisMethod)
{
    KRATOS_TRY;
    auto integration_points = rGeom.IntegrationPoints(ThisMethod);
    for (std::size_t i = 0; i < integration_points.size(); ++i)
        for (std::size_t j = 0; j < rGeom.PointsNumber(); ++j)
            TestShapeFunctionsLocalGradient(
                rGeom, integration_points[i],
                rGeom.ShapeFunctionLocalGradient(i, j, ThisMethod));
    KRATOS_CATCH("");
}

void TestShapeFunctionsLocalGradients_Point(Geometry<Node<3>> const& rGeom)
{
    KRATOS_TRY;
    Matrix local_gradient;
    auto integration_points = rGeom.IntegrationPoints();
    for (auto point : integration_points)
        TestShapeFunctionsLocalGradient(
            rGeom, point, rGeom.ShapeFunctionsLocalGradients(local_gradient, point));
    KRATOS_CATCH("");
}

void TestAllShapeFunctionsLocalGradients(Geometry<Node<3>> const& rGeom)
{
    KRATOS_TRY;
    TestShapeFunctionsLocalGradients(rGeom);
    TestShapeFunctionsLocalGradients(rGeom, GeometryData::GI_GAUSS_1);
    TestShapeFunctionsLocalGradients(rGeom, GeometryData::GI_GAUSS_2);
    TestShapeFunctionsLocalGradients_IntegrationPointIndex(rGeom);
    TestShapeFunctionsLocalGradients_IntegrationPointIndex(
        rGeom, rGeom.GetDefaultIntegrationMethod());
    TestShapeFunctionsLocalGradients_IntegrationPointIndex_ShapeFunctionIndex(
        rGeom, rGeom.GetDefaultIntegrationMethod());
    TestShapeFunctionsLocalGradients_Point(rGeom);
    KRATOS_CATCH("");
}

namespace
{
double FiniteDifference4(std::function<double(double)> f, double delta=1e-3);
}

void TestShapeFunctionsLocalGradient(Geometry<Node<3>> const& rGeom,
                                     Geometry<Node<3>>::IntegrationPointType Point,
                                     Matrix const& rLocalGradient)
{
    KRATOS_TRY;
    KRATOS_CHECK(rLocalGradient.size1() == rGeom.PointsNumber());
    KRATOS_CHECK(rLocalGradient.size2() == rGeom.LocalSpaceDimension());
    for (unsigned i = 0; i < rGeom.PointsNumber(); ++i)
        for (unsigned j = 0; j < rGeom.LocalSpaceDimension(); ++j)
        {
            auto f = [&](double delta) {
                array_1d<double, 3> point = Point;
                point[j] += delta;
                return rGeom.ShapeFunctionValue(i, point);
            };
            KRATOS_CHECK_NEAR(FiniteDifference4(f), rLocalGradient(i, j), TOLERANCE);
        }
    KRATOS_CATCH("");
}

namespace {
/**
* @brief Compute finite difference for polynomial of degree k = 4.
* @details The finite difference approximation should be exact for polynomials
* up to degree k = 4. See e.g., Mech. Struct. & Mach. 21(1), 1-66 (1993).
* @param f Scalar function of perturbation to be differentiated.
* @param delta Perturbation size.
* @return Finite difference approximation.
*/
//
double FiniteDifference4(std::function<double(double)> f, double delta)
{
    return (-f(2.0 * delta) + 8.0 * f(delta) - 8.0 * f(-delta) + f(-2.0 * delta)) /
           (12.0 * delta);
}
}
}
}
