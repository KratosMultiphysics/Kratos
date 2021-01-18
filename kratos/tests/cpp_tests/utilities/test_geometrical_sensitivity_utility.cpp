//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    
//

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/point.h"
#include "geometries/geometry.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "utilities/geometrical_sensitivity_utility.h"

namespace Kratos
{
namespace Testing
{

Geometry<Point>::Pointer CreateQuadrilateral2D4N()
{
    Geometry<Point>::PointsArrayType points;
    points.push_back(Kratos::make_shared<Point>(0.0, 0.0, 0.0));
    points.push_back(Kratos::make_shared<Point>(1.2, 0.1, 0.0));
    points.push_back(Kratos::make_shared<Point>(1.3, 1.02, 0.0));
    points.push_back(Kratos::make_shared<Point>(0.1, 1.0, 0.0));

    return Geometry<Point>::Pointer(new Quadrilateral2D4<Point>(points));
}

Geometry<Point>::Pointer CreateHexahedra3D8N()
{
    Geometry<Point>::PointsArrayType points;
    points.push_back(Kratos::make_shared<Point>(0.0, 0.0, 0.0));
    points.push_back(Kratos::make_shared<Point>(1.02, 0.0, 0.0));
    points.push_back(Kratos::make_shared<Point>(1.1, 1.03, 0.07));
    points.push_back(Kratos::make_shared<Point>(0.01, 1.0, 0.0));
    points.push_back(Kratos::make_shared<Point>(0.0, 0.0, 1.05));
    points.push_back(Kratos::make_shared<Point>(1.08, 0.01, 1.0));
    points.push_back(Kratos::make_shared<Point>(1.03, 1.09, 1.0));
    points.push_back(Kratos::make_shared<Point>(0.0, 1.05, 1.04));

    return Geometry<Point>::Pointer(new Hexahedra3D8<Point>(points));
}

Geometry<Point>::Pointer CreateTriangle2D3N()
{
    Geometry<Point>::PointsArrayType points;
    points.push_back(Kratos::make_shared<Point>(0.04, 0.02, 0.0));
    points.push_back(Kratos::make_shared<Point>(1.1, 0.03, 0.0));
    points.push_back(Kratos::make_shared<Point>(1.08, 1.0, 0.0));

    return Geometry<Point>::Pointer(new Triangle2D3<Point>(points));
}

Geometry<Point>::Pointer CreateTetrahedra3D4N()
{
    Geometry<Point>::PointsArrayType points;
    points.push_back(Kratos::make_shared<Point>(0.0, 0.04, 0.0));
    points.push_back(Kratos::make_shared<Point>(1.01, 0.0, 0.2));
    points.push_back(Kratos::make_shared<Point>(1.05, 1.0, 0.0));
    points.push_back(Kratos::make_shared<Point>(0.0, 0.03, 1.09));

    return Geometry<Point>::Pointer(new Tetrahedra3D4<Point>(points));
}

void CheckDeterminantOfJacobianSensitivityByFiniteDifference(double DetJ_Deriv,
                                                             unsigned IntegrationPoint,
                                                             ShapeParameter Deriv,
                                                             Geometry<Point>& rGeom,
                                                             GeometryData::IntegrationMethod ThisMethod,
                                                             double StepSize = 1e-7,
                                                             double Tolerance = 1e-7)
{
    double detJ = rGeom.DeterminantOfJacobian(IntegrationPoint, ThisMethod);
    rGeom[Deriv.NodeIndex].Coordinates()[Deriv.Direction] += StepSize;
    double detJ_perturbed = rGeom.DeterminantOfJacobian(IntegrationPoint, ThisMethod);
    rGeom[Deriv.NodeIndex].Coordinates()[Deriv.Direction] -= StepSize;
    double finite_difference_detJ_deriv = (detJ_perturbed - detJ) / StepSize;
    KRATOS_CHECK_NEAR(DetJ_Deriv, finite_difference_detJ_deriv, Tolerance);
}

void CheckShapeFunctionsGradientSensitivityByFiniteDifference(Matrix& DN_DX_Deriv,
                                                              unsigned IntegrationPoint,
                                                              ShapeParameter Deriv,
                                                              Geometry<Point>& rGeom,
                                                              GeometryData::IntegrationMethod ThisMethod,
                                                              double StepSize = 1e-7,
                                                              double Tolerance = 1e-7)
{
    Geometry<Point>::ShapeFunctionsGradientsType DN_DX;
    rGeom.ShapeFunctionsIntegrationPointsGradients(DN_DX, ThisMethod);
    const Matrix& rDN_DX = DN_DX[IntegrationPoint];

    rGeom[Deriv.NodeIndex].Coordinates()[Deriv.Direction] += StepSize;

    Geometry<Point>::ShapeFunctionsGradientsType DN_DX_perturbed;
    rGeom.ShapeFunctionsIntegrationPointsGradients(DN_DX_perturbed, ThisMethod);
    const Matrix& rDN_DX_perturbed = DN_DX_perturbed[IntegrationPoint];

    rGeom[Deriv.NodeIndex].Coordinates()[Deriv.Direction] -= StepSize;

    // Perform some checks to avoid false positives.
    KRATOS_CHECK(DN_DX_Deriv.size1() != 0);
    KRATOS_CHECK(DN_DX_Deriv.size2() != 0);
    for (unsigned i = 0; i < DN_DX_Deriv.size1(); ++i)
        for (unsigned j = 0; j < DN_DX_Deriv.size2(); ++j)
            {
                const double finite_difference_ij = (rDN_DX_perturbed(i,j) - rDN_DX(i,j)) / StepSize;
                KRATOS_CHECK_NEAR(DN_DX_Deriv(i, j), finite_difference_ij, Tolerance);
            }
}

void TestThisGeometry(Geometry<Point>& rGeom,
                      GeometryData::IntegrationMethod ThisMethod,
                      double StepSize = 1e-7,
                      double Tolerance = 1e-7)
{

    Geometry<Point>::JacobiansType J;
    rGeom.Jacobian(J, ThisMethod);

    Geometry<Point>::ShapeFunctionsGradientsType DN_De;
    DN_De = rGeom.ShapeFunctionsLocalGradients(ThisMethod);

    GeometricalSensitivityUtility::ShapeFunctionsGradientType DN_DX_deriv;

    for (unsigned g = 0; g < rGeom.IntegrationPointsNumber(ThisMethod); ++g)
    {
        const Matrix& rJ = J[g];
        const Matrix& rDN_De = DN_De[g];
        GeometricalSensitivityUtility geom_sensitivity(rJ, rDN_De);

        auto s = ShapeParameter::Sequence(rGeom.PointsNumber(),
                                          rGeom.WorkingSpaceDimension());
        while (s)
        {
            double detJ_deriv;
            const auto& deriv = s.CurrentValue();
            geom_sensitivity.CalculateSensitivity(deriv, detJ_deriv, DN_DX_deriv);
            CheckDeterminantOfJacobianSensitivityByFiniteDifference(
                detJ_deriv, g, deriv, rGeom, ThisMethod, StepSize, Tolerance);
            CheckShapeFunctionsGradientSensitivityByFiniteDifference(
                DN_DX_deriv, g, deriv, rGeom, ThisMethod, StepSize, Tolerance);
            ++s;
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(GeometricalSensitivityUtility_Quadrilateral2D4N_GAUSS_2, KratosSensitivityTestSuite)
{
    Geometry<Point>::Pointer p_geom = CreateQuadrilateral2D4N();
    Geometry<Point>& r_geom = *p_geom;
    TestThisGeometry(r_geom, GeometryData::GI_GAUSS_2);
}

KRATOS_TEST_CASE_IN_SUITE(GeometricalSensitivityUtility_Quadrilateral2D4N_GAUSS_1, KratosSensitivityTestSuite)
{
    Geometry<Point>::Pointer p_geom = CreateQuadrilateral2D4N();
    Geometry<Point>& r_geom = *p_geom;
    TestThisGeometry(r_geom, GeometryData::GI_GAUSS_1);
}

KRATOS_TEST_CASE_IN_SUITE(GeometricalSensitivityUtility_Hexahedra3D8N_GAUSS_2, KratosSensitivityTestSuite)
{
    Geometry<Point>::Pointer p_geom = CreateHexahedra3D8N();
    Geometry<Point>& r_geom = *p_geom;
    TestThisGeometry(r_geom, GeometryData::GI_GAUSS_2);
}

KRATOS_TEST_CASE_IN_SUITE(GeometricalSensitivityUtility_Triangle2D3N_GAUSS_1, KratosSensitivityTestSuite)
{
    Geometry<Point>::Pointer p_geom = CreateTriangle2D3N();
    Geometry<Point>& r_geom = *p_geom;
    TestThisGeometry(r_geom, GeometryData::GI_GAUSS_1, 1e-8);
}

KRATOS_TEST_CASE_IN_SUITE(GeometricalSensitivityUtility_Tetrahedra3D4N_GAUSS_1, KratosSensitivityTestSuite)
{
    Geometry<Point>::Pointer p_geom = CreateTetrahedra3D4N();
    Geometry<Point>& r_geom = *p_geom;
    TestThisGeometry(r_geom, GeometryData::GI_GAUSS_1, 1e-8);
}

} // namespace Testing
} // namespace Kratos.