// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include <algorithm>

#include "math_utilities.h"
#include "utilities/math_utils.h"

namespace Kratos
{

std::vector<double> GeoMechanicsMathUtilities::CalculateDeterminants(const std::vector<Matrix>& rMatrices)
{
    std::vector<double> result(rMatrices.size());
    std::transform(rMatrices.cbegin(), rMatrices.cend(), result.begin(),
                   [](const auto& rMatrix) { return MathUtils<>::Det(rMatrix); });

    return result;
}

Matrix GeoMechanicsMathUtilities::VectorToDiagonalMatrix(const Vector& rVector)
{
    Matrix result = ZeroMatrix(rVector.size(), rVector.size());
    for (std::size_t i = 0; i < rVector.size(); ++i) {
        result(i, i) = rVector(i);
    }
    return result;
}

Vector GeoMechanicsMathUtilities::DiagonalMatrixToVector(const Matrix& rMatrix)
{
    KRATOS_ERROR_IF(rMatrix.size1() != rMatrix.size2())
        << "Attempting to convert diagonal matrix to vector, but the matrix is not square\n";
    bool found_non_diagonal = false;
    for (std::size_t i = 0; i < rMatrix.size1(); ++i) {
        for (std::size_t j = 0; j < rMatrix.size2(); ++j) {
            if (rMatrix(i, j) != 0.0 && i != j) {
                KRATOS_WARNING("DiagonalMatrixToVector")
                    << "The matrix is not diagonal. A non-zero off-diagonal component is found at ("
                    << i << ", " << j << ") with value " << rMatrix(i, j) << "." << std::endl;
                found_non_diagonal = true;
                break;
            }
        }
        if (found_non_diagonal) break;
    }

    Vector result = ZeroVector(rMatrix.size1());
    for (std::size_t i = 0; i < rMatrix.size1(); ++i) {
        result(i) = rMatrix(i, i);
    }
    return result;
}

Matrix GeoMechanicsMathUtilities::RotateSecondOrderTensor(const Matrix& rTensor, const Matrix& rRotationMatrix)
{
    return prod(rRotationMatrix, Matrix{prod(rTensor, trans(rRotationMatrix))});
}

} // namespace Kratos