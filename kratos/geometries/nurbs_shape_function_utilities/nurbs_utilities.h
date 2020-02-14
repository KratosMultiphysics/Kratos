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

// System includes

// External includes

// Project includes

namespace Kratos {

///@name Kratos Classes
///@{

/// Utility functions for NURBS computation
/*
* Provides universal geometrical utiltity functions for the computation of
* curve and surface NURBS/ B-Spline shape functions.
 */
namespace NurbsUtilities
{
    ///@name Type Definitions
    ///@{

    typedef typename std::size_t IndexType;
    typedef typename std::size_t SizeType;

    ///@}
    ///@name Static Operations
    ///@{
    /*
    * @brief the index of the upper limit of the span in which the ParameterT lays.
    * From Piegl and Tiller, The NURBS Book, Algorithm A2.1
    */
    inline IndexType GetUpperSpan(
        const SizeType PolynomialDegree,
        const Vector& rKnots,
        const double ParameterT)
    {
        const auto span = std::upper_bound(std::begin(rKnots) + PolynomialDegree,
            std::end(rKnots) - PolynomialDegree, ParameterT) - std::begin(rKnots) - 1;
        return span;
    }

    /*
    * @brief the index of the lower limit of the span in which the ParameterT lays.
    * From Piegl and Tiller, The NURBS Book, Algorithm A2.1
    */
    inline IndexType GetLowerSpan(
        const SizeType PolynomialDegree,
        const Vector& rKnots,
        const double ParameterT)
    {
        const auto span = std::lower_bound(std::begin(rKnots) + PolynomialDegree,
            std::end(rKnots) - PolynomialDegree, ParameterT) - std::begin(rKnots) - 1;
        return span;
    }

    /*
    * @brief Computes the degree of a nurbs/ b-spline shape by:
    * @param NumberOfKnots and
    * @param NumberOfControlPoints
    */
    inline SizeType GetPolynomialDegree(const SizeType NumberOfKnots, const SizeType NumberOfControlPoints)
    {
        return NumberOfKnots - NumberOfControlPoints + 1;
    }

    /*
    * @brief Computes the number of knots of a nurbs/ b-spline shape by:
    * @param PolynomialDegree and
    * @param NumberOfControlPoints
    */
    inline SizeType GetNumberOfKnots(const SizeType PolynomialDegree, const SizeType NumberOfControlPoints)
    {
        return NumberOfControlPoints + PolynomialDegree - 1;
    }

    /*
    * @brief Computes the number of control points of a nurbs/ b-spline shape by:
    * @param PolynomialDegree and
    * @param NumberOfKnots
    */
    inline SizeType GetNumberOfControlPoints(const SizeType PolynomialDegree, const SizeType NumberOfKnots)
    {
        return NumberOfKnots - PolynomialDegree + 1;
    }

    /*
    * @brief Computes the number of spans of a nurbs/ b-spline shape by:
    * @param PolynomialDegree and
    * @param NumberOfKnots
    */
    inline SizeType GetNumberOfSpans(const SizeType PolynomialDegree, const SizeType NumberOfKnots)
    {
        return NumberOfKnots - 2 * PolynomialDegree + 1;
    }

    /*
    * @brief Computes the binomial coefficient for (N || K).
    */
    static constexpr inline SizeType GetBinomCoefficient(const SizeType N, const SizeType K) noexcept
    {
        return
            (K > N) ? 0 :  // out of range
            (K == 0 || K == N) ? 1 :  // edge
            (K == 1 || K == N - 1) ? N :  // first
            (K + K < N) ?                // recursive:
            (GetBinomCoefficient(N - 1, K - 1) * N) / K :  //   path to K = 1     faster
            (GetBinomCoefficient(N - 1, K) * N) / (N - K);  //   path to K = n - 1 faster
    }


    /*
    * @brief Computes a vector index from two matrix indicies.
    * @return index within vector
    */
    static constexpr inline IndexType GetVectorIndexFromMatrixIndices(
        const SizeType NumberPerRow, const SizeType NumberPerColumn,
        const IndexType RowIndex, const IndexType ColumnIndex) noexcept
    {
        return ColumnIndex * NumberPerRow + RowIndex;
    }

    /*
    * @brief Computes two matrix indices from vector index.
    * @return indices within Matrix
    */
    static inline std::pair<IndexType, IndexType> GetMatrixIndicesFromVectorIndex(
        const SizeType NumberPerRow,
        const SizeType NumberPerColumn,
        const IndexType Index) noexcept
    {
        const IndexType row = Index % NumberPerRow;
        const IndexType col = Index / NumberPerRow;

        return std::make_pair(row, col);
    }
    ///@}
}; // class NurbsUtility
///@}
} // namespace Kratos

#endif // KRATOS_NURBS_UTILITY_H_INCLUDED defined
