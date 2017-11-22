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
#include <vector>

// Project includes
#include "testing/testing.h"
#include "geometries/point.h"
#include "geometries/geometry.h"
#include "geometries/quadrilateral_2d_4.h"

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

KRATOS_TEST_CASE_IN_SUITE(GeometricalSensitivityUtility_Constructor, KratosSensitivityTestSuite)
{
    const std::size_t n = 2;
    double val = 1.0;
    GeometricalSensitivityUtility::JacobianType J(n, n);
    for (unsigned i = 0; i < n; ++i)
        for (unsigned j = 0; j < n; ++j)
        {
            J(i, j) = val;
            val += 1.0;
        }
    GeometricalSensitivityUtility::ShapeFunctionsLocalGradientType DN_De(n, n);
    GeometricalSensitivityUtility geom_sensitivity(J, DN_De);
}

KRATOS_TEST_CASE_IN_SUITE(GeometricalSensitivityUtility_Quadrilateral2D4N, KratosSensitivityTestSuite)
{
    Geometry<Point>::Pointer p_geom = CreateQuadrilateral2D4N();
    Geometry<Point>& r_geom = *p_geom;
    Geometry<Point>::JacobiansType J;
    Geometry<Point>::ShapeFunctionsGradientsType DN_De;
    r_geom.Jacobian(J);
    DN_De = r_geom.ShapeFunctionsLocalGradients();
    GeometricalSensitivityUtility::ShapeFunctionsGradientType DN_DX_deriv;
    for (unsigned g = 0; g < J.size(); ++g)
    {
        GeometricalSensitivityUtility geom_sensitivity(J(g), DN_De(g));
        std::vector<std::vector<double>> detJ_deriv(r_geom.size());
        for (unsigned i_node = 0; i_node < r_geom.size(); ++i_node)
            detJ_deriv[i_node].resize((*p_geom)[i_node].Dimension());
        for (unsigned i_node = 0; i_node < p_geom->size(); ++i_node)
            for (unsigned i_coord = 0; i_coord < r_geom[i_node].Dimension(); i_coord++)
            {
                geom_sensitivity.CalculateSensitivity(
                    i_node, i_coord, detJ_deriv[i_node][i_coord], DN_DX_deriv);
                std::cout << '(' << i_node << ',' << i_coord
                          << ") : " << detJ_deriv[i_node][i_coord] << std::endl;
            }
        // CheckFiniteDifferenceSensitivity(p_geom, detJ_deriv);
    }
}

} // namespace Testing
} // namespace Kratos.