//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: AdjointFluidApplication/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//

// External includes
//#include <iomanip>

// Project includes
#include "testing/testing.h"
#include "geometries/point.h"
#include "geometries/geometry.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/hexahedra_3d_8.h"

// Application includes
#include "custom_utilities/geometrical_sensitivity_utility.h"

namespace Kratos
{
namespace Testing
{

Geometry<Point>::Pointer CreateQuadrilateral2D4N()
{
    Geometry<Point>::PointsArrayType points;
    points.push_back(boost::make_shared<Point>(0.0, 0.0, 0.0));
    points.push_back(boost::make_shared<Point>(1.0, 0.0, 0.0));
    points.push_back(boost::make_shared<Point>(1.0, 1.0, 0.0));
    points.push_back(boost::make_shared<Point>(0.0, 1.0, 0.0));

    return Geometry<Point>::Pointer(new Quadrilateral2D4<Point>(points));
}

Geometry<Point>::Pointer CreateHexahedra3D8N()
{
    Geometry<Point>::PointsArrayType points;
    points.push_back(boost::make_shared<Point>(0.0, 0.0, 0.0));
    points.push_back(boost::make_shared<Point>(1.0, 0.0, 0.0));
    points.push_back(boost::make_shared<Point>(1.0, 1.0, 0.0));
    points.push_back(boost::make_shared<Point>(0.0, 1.0, 0.0));
    points.push_back(boost::make_shared<Point>(0.0, 0.0, 1.0));
    points.push_back(boost::make_shared<Point>(1.0, 0.0, 1.0));
    points.push_back(boost::make_shared<Point>(1.0, 1.0, 1.0));
    points.push_back(boost::make_shared<Point>(0.0, 1.0, 1.0));

    return Geometry<Point>::Pointer(new Hexahedra3D8<Point>(points));
}

void CheckDeterminantOfJacobianSensitivityByFiniteDifference(double DetJ_Deriv,
                                                             unsigned IntegrationPoint,
                                                             unsigned iNode,
                                                             unsigned iCoord,
                                                             Geometry<Point>& rGeom,
                                                             GeometryData::IntegrationMethod ThisMethod,
                                                             double StepSize = 1e-7,
                                                             double Tolerance = 1e-7)
{
    double detJ = rGeom.DeterminantOfJacobian(IntegrationPoint, ThisMethod);
    rGeom[iNode].Coordinates()[iCoord] += StepSize;
    double detJ_perturbed = rGeom.DeterminantOfJacobian(IntegrationPoint, ThisMethod);
    rGeom[iNode].Coordinates()[iCoord] -= StepSize;
    double finite_difference_detJ_deriv = (detJ_perturbed - detJ) / StepSize;
    KRATOS_CHECK_NEAR(DetJ_Deriv, finite_difference_detJ_deriv, Tolerance);
}

void CheckShapeFunctionsGradientSensitivityByFiniteDifference(Matrix& DN_DX_Deriv,
                                                              unsigned IntegrationPoint,
                                                              unsigned iNode,
                                                              unsigned iCoord,
                                                              Geometry<Point>& rGeom,
                                                              GeometryData::IntegrationMethod ThisMethod,
                                                              double StepSize = 1e-7,
                                                              double Tolerance = 1e-7)
{
    Geometry<Point>::ShapeFunctionsGradientsType DN_DX;
    rGeom.ShapeFunctionsIntegrationPointsGradients(DN_DX, ThisMethod);
    const Matrix& rDN_DX = DN_DX[IntegrationPoint];

    rGeom[iNode].Coordinates()[iCoord] += StepSize;

    Geometry<Point>::ShapeFunctionsGradientsType DN_DX_perturbed;
    rGeom.ShapeFunctionsIntegrationPointsGradients(DN_DX_perturbed, ThisMethod);
    const Matrix& rDN_DX_perturbed = DN_DX_perturbed[IntegrationPoint];

    rGeom[iNode].Coordinates()[iCoord] -= StepSize;

    for (unsigned i = 0; i < DN_DX_Deriv.size1(); ++i)
        for (unsigned j = 0; j < DN_DX_Deriv.size2(); ++j)
            {
                const double finite_difference_ij = (rDN_DX_perturbed(i,j) - rDN_DX(i,j)) / StepSize;
                KRATOS_CHECK_NEAR(DN_DX_Deriv(i,j), finite_difference_ij, Tolerance);
            }
}

KRATOS_TEST_CASE_IN_SUITE(GeometricalSensitivityUtility_Quadrilateral2D4N_GAUSS_2, KratosSensitivityTestSuite)
{
    Geometry<Point>::Pointer p_geom = CreateQuadrilateral2D4N();
    Geometry<Point>& r_geom = *p_geom;

    // Geometry<Point>::IntegrationPointsArrayType gauss_points =
    //     r_geom.IntegrationPoints(GeometryData::GI_GAUSS_2);
    // std::cout << "Shape function values : ";
    // for (unsigned i = 0; i < gauss_points.size(); ++i)
    //     std::cout << gauss_points[i] << std::endl;

    Geometry<Point>::JacobiansType J;
    r_geom.Jacobian(J, GeometryData::GI_GAUSS_2);

    Geometry<Point>::ShapeFunctionsGradientsType DN_De;
    DN_De = r_geom.ShapeFunctionsLocalGradients(GeometryData::GI_GAUSS_2);

    GeometricalSensitivityUtility::ShapeFunctionsGradientType DN_DX_deriv;

    for (unsigned g = 0; g < r_geom.IntegrationPointsNumber(GeometryData::GI_GAUSS_2); ++g)
    {
        const Matrix& rJ = J[g];
        const Matrix& rDN_De = DN_De[g];
        GeometricalSensitivityUtility geom_sensitivity(rJ, rDN_De);

        for (unsigned i_node = 0; i_node < p_geom->size(); ++i_node)
            for (unsigned i_coord = 0; i_coord < r_geom[i_node].Dimension(); i_coord++)
            {
                double detJ_deriv;
                geom_sensitivity.CalculateSensitivity(i_node, i_coord, detJ_deriv, DN_DX_deriv);
                CheckDeterminantOfJacobianSensitivityByFiniteDifference(
                    detJ_deriv, g, i_node, i_coord, r_geom, GeometryData::GI_GAUSS_2);
                CheckShapeFunctionsGradientSensitivityByFiniteDifference(
                    DN_DX_deriv, g, i_node, i_coord, r_geom, GeometryData::GI_GAUSS_2);
            }
    }
}

KRATOS_TEST_CASE_IN_SUITE(GeometricalSensitivityUtility_Hexahedra3D8N_GAUSS_2, KratosSensitivityTestSuite)
{
    Geometry<Point>::Pointer p_geom = CreateHexahedra3D8N();
    Geometry<Point>& r_geom = *p_geom;

    Geometry<Point>::JacobiansType J;
    r_geom.Jacobian(J, GeometryData::GI_GAUSS_2);

    Geometry<Point>::ShapeFunctionsGradientsType DN_De;
    DN_De = r_geom.ShapeFunctionsLocalGradients(GeometryData::GI_GAUSS_2);

    GeometricalSensitivityUtility::ShapeFunctionsGradientType DN_DX_deriv;

    for (unsigned g = 0; g < r_geom.IntegrationPointsNumber(GeometryData::GI_GAUSS_2); ++g)
    {
        const Matrix& rJ = J[g];
        const Matrix& rDN_De = DN_De[g];
        GeometricalSensitivityUtility geom_sensitivity(rJ, rDN_De);

        for (unsigned i_node = 0; i_node < p_geom->size(); ++i_node)
            for (unsigned i_coord = 0; i_coord < r_geom[i_node].Dimension(); i_coord++)
            {
                double detJ_deriv;
                geom_sensitivity.CalculateSensitivity(i_node, i_coord, detJ_deriv, DN_DX_deriv);
                CheckDeterminantOfJacobianSensitivityByFiniteDifference(
                    detJ_deriv, g, i_node, i_coord, r_geom, GeometryData::GI_GAUSS_2);
                CheckShapeFunctionsGradientSensitivityByFiniteDifference(
                    DN_DX_deriv, g, i_node, i_coord, r_geom, GeometryData::GI_GAUSS_2);
            }
    }
}

} // namespace Testing
} // namespace Kratos.