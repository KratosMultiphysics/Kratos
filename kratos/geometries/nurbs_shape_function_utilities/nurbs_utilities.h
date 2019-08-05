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

#if !defined(KRATOS_NURBS_UTILITY_H_INCLUDED )
#define  KRATOS_NURBS_UTILITY_H_INCLUDED

#include <algorithm>

namespace Kratos {

///@name Kratos Classes
///@{
/// 
/** Provides universal geometrical utiltity functions for the computation of
    curve and surface NURBS/ B-Spline shape functions.
 */
class NurbsUtilities
{
public:
    ///@name Static Operations
    ///@{

    static int GetUpperSpan(
        const int PolynomialDegree,
        const Vector& rKnots,
        const double ParameterT)
    {
        auto span = std::upper_bound(std::begin(rKnots) + PolynomialDegree,
            std::end(rKnots) - PolynomialDegree, ParameterT) - std::begin(rKnots) - 1;
        return static_cast<int>(span);
    }

    static int GetLowerSpan(
        const int PolynomialDegree,
        const Vector& rKnots,
        const double ParameterT)
    {
        auto span = std::lower_bound(std::begin(rKnots) + PolynomialDegree,
            std::end(rKnots) - PolynomialDegree, ParameterT) - std::begin(rKnots) - 1;
        return static_cast<int>(span);
    }

    /* Computes the degree of a nurbs/ b-spline shape by:
    @param NumberOfKnots and
    @param NumberOfControlPoints*/
    static int GetPolynomialDegree(const int NumberOfKnots, const int NumberOfControlPoints)
    {
        return NumberOfKnots - NumberOfControlPoints + 1;
    }

    /* Computes the number of knots of a nurbs/ b-spline shape by:
    @param PolynomialDegree and
    @param NumberOfControlPoints*/
    static int GetNumberOfKnots(const int PolynomialDegree, const int NumberOfControlPoints)
    {
        return NumberOfControlPoints + PolynomialDegree - 1;
    }

    /* Computes the number of control points of a nurbs/ b-spline shape by:
    @param PolynomialDegree and
    @param NumberOfKnots*/
    static int GetNumberOfControlPoints(const int PolynomialDegree, const int NumberOfKnots)
    {
        return NumberOfKnots - PolynomialDegree + 1;
    }

    /* Computes the number of spans of a nurbs/ b-spline shape by:
    @param PolynomialDegree and
    @param NumberOfKnots*/
    static int GetNumberOfSpans(const int PolynomialDegree, const int NumberOfKnots)
    {
        return NumberOfKnots - 2 * PolynomialDegree + 1;
    }

    /* Computes the binomial coefficient for (N || K).*/
    static constexpr inline int GetBinomCoefficient(const int N, const int K) noexcept
    {
        // clang-format off
        return
            (K > N) ? 0 :  // out of range
            (K == 0 || K == N) ? 1 :  // edge
            (K == 1 || K == N - 1) ? N :  // first
            (K + K < N) ?                // recursive:
            (GetBinomCoefficient(N - 1, K - 1) * N) / K :  //   path to K = 1     faster
            (GetBinomCoefficient(N - 1, K) * N) / (N - K);  //   path to K = n - 1 faster
        // clang-format on
    }


    /* Computes the vector index from the matrix indicies*/
    static constexpr inline int GetVectorIndexFromMatrixIndices(
        const int Rows, const int Columns,
        const int Row, const int Column) noexcept
    {
        return Row * Columns + Column;
    }

    /* Computes the matrix indices from the vector index*/
    static inline std::pair<int, int> GetMatrixIndicesFromVectorIndex(
        const int Rows,
        const int Columns,
        const int Index) noexcept
    {
        const int row = Index / Columns;
        const int col = Index % Columns;

        return std::make_pair(row, col);
    }
    ///@}
}; // class NurbsUtility
///@}
} // namespace Kratos

#endif // KRATOS_NURBS_UTILITY_H_INCLUDED defined 