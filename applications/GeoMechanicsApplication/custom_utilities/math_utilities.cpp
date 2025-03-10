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

Matrix GeoMechanicsMathUtilities::RotateTensor(const Matrix& rTensor, const Matrix& rRotationMatrix)
{
    Matrix temp = prod(rTensor, trans(rRotationMatrix));
    return prod(rRotationMatrix, temp);
}

} // namespace Kratos