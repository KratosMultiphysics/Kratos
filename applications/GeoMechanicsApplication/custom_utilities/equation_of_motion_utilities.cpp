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

Vector GeoEquationOfMotionUtilities::CalculateInternalForceVector(const std::vector<Matrix>& rBs,
                                                                  const std::vector<Vector>& rStressVectors,
                                                                  const std::vector<double>& rIntegrationCoefficients)
{
    KRATOS_ERROR_IF((rBs.size() != rStressVectors.size()) || (rBs.size() != rIntegrationCoefficients.size()))
        << "Cannot calculate the internal force vector: input vectors have different sizes\n";
    KRATOS_ERROR_IF(rBs.empty())
        << "Cannot calculate the internal force vector: input vectors are empty\n";
    auto has_inconsistent_sizes = [number_of_rows    = rBs.front().size1(),
                                   number_of_columns = rBs.front().size2()](const auto& rMatrix) {
        return (rMatrix.size1() != number_of_rows) || (rMatrix.size2() != number_of_columns);
    };
    KRATOS_ERROR_IF(std::any_of(rBs.begin() + 1, rBs.end(), has_inconsistent_sizes))
        << "Cannot calculate the internal force vector: B-matrices have different sizes";

    auto result = Vector{ZeroVector{rBs.front().size2()}};
    for (auto i = std::size_t{0}; i < rBs.size(); ++i) {
        result += prod(trans(rBs[i]), rStressVectors[i]) * rIntegrationCoefficients[i];
    }

    return result;
}

} /* namespace Kratos.*/
