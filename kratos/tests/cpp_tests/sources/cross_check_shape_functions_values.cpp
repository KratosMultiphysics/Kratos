#include "tests/cpp_tests/geometries/cross_check_shape_functions_values.h"

//#include "containers/array_1d.h"
#include "testing/testing.h"
#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

namespace {
void Check_ShapeFunctionsValues1(Geometry<Node<3>> const& rGeom);
void Check_ShapeFunctionsValues2(Geometry<Node<3>> const& rGeom);
void Check_ShapeFunctionsValues3(Geometry<Node<3>> const& rGeom);
void Check_ShapeFunctionsValues4(Geometry<Node<3>> const& rGeom);
void Check_ShapeFunctionsValues5(Geometry<Node<3>> const& rGeom);
void Check(double N,
           Geometry<Node<3>> const& rGeom,
           Geometry<Node<3>>::IndexType ShapeFunctionIndex,
           array_1d<double, 3> const& rCoord);
}

void CrossCheckShapeFunctionsValues(Geometry<Node<3>> const& rGeom)
{
    KRATOS_TRY;
    Check_ShapeFunctionsValues1(rGeom);
    Check_ShapeFunctionsValues2(rGeom);
    Check_ShapeFunctionsValues3(rGeom);
    Check_ShapeFunctionsValues4(rGeom);
    Check_ShapeFunctionsValues5(rGeom);
    KRATOS_CATCH("");
}

namespace {
void Check_ShapeFunctionsValues1(Geometry<Node<3>> const& rGeom)
{
    KRATOS_TRY;
    const Matrix& r_shape_functions_values = rGeom.ShapeFunctionsValues();
    const auto& r_integration_points = rGeom.IntegrationPoints();
    for (std::size_t g = 0; g < rGeom.IntegrationPointsNumber(); ++g)
        for (std::size_t i = 0; i < rGeom.PointsNumber(); ++i)
            Check(r_shape_functions_values(g, i), rGeom, i, r_integration_points[g]);
    KRATOS_CATCH("");
}

void Check_ShapeFunctionsValues2(Geometry<Node<3>> const& rGeom)
{
    KRATOS_TRY;
    const auto& r_integration_points = rGeom.IntegrationPoints();
    Vector shape_functions_values(rGeom.PointsNumber());
    for (std::size_t g = 0; g < rGeom.IntegrationPointsNumber(); ++g)
    {
        rGeom.ShapeFunctionsValues(shape_functions_values, r_integration_points[g]);
        for (std::size_t i = 0; i < rGeom.PointsNumber(); ++i)
            Check(shape_functions_values[i], rGeom, i, r_integration_points[g]);
    }
    KRATOS_CATCH("");
}

void Check_ShapeFunctionsValues3(Geometry<Node<3>> const& rGeom)
{
    KRATOS_TRY;
    const Matrix& r_shape_functions_values =
        rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);
    const auto& r_integration_points = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_1);
    for (std::size_t g = 0; g < rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_1); ++g)
        for (std::size_t i = 0; i < rGeom.PointsNumber(); ++i)
            Check(r_shape_functions_values(g, i), rGeom, i, r_integration_points[g]);
    KRATOS_CATCH("");
}

void Check_ShapeFunctionsValues4(Geometry<Node<3>> const& rGeom)
{
    KRATOS_TRY;
    const auto& r_integration_points = rGeom.IntegrationPoints();
    for (std::size_t g = 0; g < rGeom.IntegrationPointsNumber(); ++g)
        for (std::size_t i = 0; i < rGeom.PointsNumber(); ++i)
            Check(rGeom.ShapeFunctionValue(g, i), rGeom, i, r_integration_points[g]);
    KRATOS_CATCH("");
}

void Check_ShapeFunctionsValues5(Geometry<Node<3>> const& rGeom)
{
    KRATOS_TRY;
    const auto& r_integration_points = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_1);
    for (std::size_t g = 0; g < rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_1); ++g)
        for (std::size_t i = 0; i < rGeom.PointsNumber(); ++i)
            Check(rGeom.ShapeFunctionValue(g, i, GeometryData::GI_GAUSS_1), rGeom, i, r_integration_points[g]);
    KRATOS_CATCH("");
}

void Check(double N,
           Geometry<Node<3>> const& rGeom,
           Geometry<Node<3>>::IndexType ShapeFunctionIndex,
           array_1d<double, 3> const& rCoord)
{
    KRATOS_TRY;
    KRATOS_CHECK_NEAR(N, rGeom.ShapeFunctionValue(ShapeFunctionIndex, rCoord), TOLERANCE);
    KRATOS_CATCH("");
}

}
}
}
