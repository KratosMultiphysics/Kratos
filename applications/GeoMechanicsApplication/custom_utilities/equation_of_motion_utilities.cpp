// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

// Project includes

// Application includes
#include "custom_utilities/equation_of_motion_utilities.h"
#include "custom_utilities/element_utilities.hpp"
#include "utilities/geometry_utilities.h"

namespace Kratos
{

Matrix GeoEquationOfMotionUtilities::CalculateMassMatrix(std::size_t   dimension,
                                                         std::size_t   number_U_nodes,
                                                         std::size_t   NumberIntegrationPoints,
                                                         const Matrix& Nu_container,
                                                         const std::vector<double>& rSolidDensities,
                                                         const std::vector<double>& rIntegrationCoefficients)
{
    const std::size_t block_element_size = number_U_nodes * dimension;
    Matrix            Nu                 = ZeroMatrix(dimension, block_element_size);
    Matrix            mass_matrix        = ZeroMatrix(block_element_size, block_element_size);

    for (unsigned int g_point = 0; g_point < NumberIntegrationPoints; ++g_point) {
        GeoElementUtilities::CalculateNuMatrix(dimension, number_U_nodes, Nu, Nu_container, g_point);

        mass_matrix += rSolidDensities[g_point] * prod(trans(Nu), Nu) * rIntegrationCoefficients[g_point];
    }
    return mass_matrix;
}

Vector GeoEquationOfMotionUtilities::CalculateDetJsInitialConfiguration(const Geometry<Node>& rGeom,
                                                                        const GeometryData::IntegrationMethod IntegrationMethod)
{
    const Geometry<Node>::IntegrationPointsArrayType& integration_points =
        rGeom.IntegrationPoints(IntegrationMethod);
    const std::size_t number_G_points = integration_points.size();

    Vector det_Js_initial_configuration(number_G_points);
    Matrix J0;
    Matrix inv_J0;
    for (unsigned int g_point = 0; g_point < number_G_points; ++g_point) {
        GeometryUtils::JacobianOnInitialConfiguration(rGeom, integration_points[g_point], J0);
        MathUtils<double>::InvertMatrix(J0, inv_J0, det_Js_initial_configuration(g_point));
    }
    return det_Js_initial_configuration;
}

Matrix GeoEquationOfMotionUtilities::CalculateDampingMatrix(double        RayleighAlpha,
                                                            double        RayleighBeta,
                                                            const Matrix& rMassMatrix,
                                                            const Matrix& rStiffnessMatrix)
{
    return RayleighAlpha * rMassMatrix + RayleighBeta * rStiffnessMatrix;
}

Matrix GeoEquationOfMotionUtilities::CalculateStiffnessMatrixGPoint(const Matrix& rB,
                                                                    const Matrix& rConstitutiveMatrix,
                                                                    double IntegrationCoefficient)
{
    return prod(trans(rB), Matrix(prod(rConstitutiveMatrix, rB))) * IntegrationCoefficient;
}

Matrix GeoEquationOfMotionUtilities::CalculateStiffnessMatrix(const std::vector<Matrix>& rBs,
                                                              const std::vector<Matrix>& rConstitutiveMatrices,
                                                              const std::vector<double>& rIntegrationCoefficients)
{
    Matrix result = ZeroMatrix(rBs[0].size2(), rBs[0].size2());
    for (unsigned int GPoint = 0; GPoint < rBs.size(); ++GPoint) {
        result += CalculateStiffnessMatrixGPoint(rBs[GPoint], rConstitutiveMatrices[GPoint],
                                                 rIntegrationCoefficients[GPoint]);
    }
    return result;
}

} /* namespace Kratos.*/
