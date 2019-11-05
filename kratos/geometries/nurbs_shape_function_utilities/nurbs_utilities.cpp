//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Thomas Oberbichler
//
//  Ported from the ANurbs library (https://github.com/oberbichler/ANurbs)
//

// System includes

// External includes

// Project includes
#include "geometries/nurbs_shape_function_utilities/nurbs_utilities.h"

namespace Kratos {

namespace NurbsUtilities
{

std::size_t GetUpperSpan(
    const std::size_t PolynomialDegree,
    const Vector& rKnots,
    const double ParameterT
    )
{
    const auto span = std::upper_bound(std::begin(rKnots) + PolynomialDegree,
        std::end(rKnots) - PolynomialDegree, ParameterT) - std::begin(rKnots) - 1;
    return span;
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t GetLowerSpan(
    const std::size_t PolynomialDegree,
    const Vector& rKnots,
    const double ParameterT
    )
{
    const auto span = std::lower_bound(std::begin(rKnots) + PolynomialDegree,
        std::end(rKnots) - PolynomialDegree, ParameterT) - std::begin(rKnots) - 1;
    return span;
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t GetPolynomialDegree(
    const std::size_t NumberOfKnots, 
    const std::size_t NumberOfControlPoints
    )
{
    return NumberOfKnots - NumberOfControlPoints + 1;
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t GetNumberOfKnots(
    const std::size_t PolynomialDegree, 
    const std::size_t NumberOfControlPoints
    )
{
    return NumberOfControlPoints + PolynomialDegree - 1;
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t GetNumberOfControlPoints(
    const std::size_t PolynomialDegree, 
    const std::size_t NumberOfKnots
    )
{
    return NumberOfKnots - PolynomialDegree + 1;
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t GetNumberOfSpans(
    const std::size_t PolynomialDegree, 
    const std::size_t NumberOfKnots
    )
{
    return NumberOfKnots - 2 * PolynomialDegree + 1;
}

/***********************************************************************************/
/***********************************************************************************/

std::pair<std::size_t, std::size_t> GetMatrixIndicesFromVectorIndex(
    const std::size_t NumberPerRow,
    const std::size_t NumberPerColumn,
    const std::size_t Index
    ) noexcept
{
    const IndexType row = Index % NumberPerRow;
    const IndexType col = Index / NumberPerRow;

    return std::make_pair(row, col);
}
}; // class NurbsUtility
} // namespace Kratos
