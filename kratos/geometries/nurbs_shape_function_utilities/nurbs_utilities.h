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
#include <utility>
#include <algorithm>
#include <iterator>

// External includes

// Project includes
#include "includes/define.h"
#include "spaces/ublas_space.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @namespace NurbsUtilities
 * @ingroup KratosCore
 * @brief Utility functions for NURBS computation
 * @details Provides universal geometrical utiltity functions for the computation of curve and surface NURBS/ B-Spline shape functions.
 * @author Thomas Oberbichler
 */
namespace NurbsUtilities
{
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    ///@}
    ///@name Static Operations
    ///@{
    
    /**
     * @brief the index of the upper limit of the span in which the ParameterT lays.
     * @note From Piegl and Tiller, The NURBS Book, Algorithm A2.1
     */
    IndexType KRATOS_API(KRATOS_CORE) GetUpperSpan(
        const SizeType PolynomialDegree,
        const Vector& rKnots,
        const double ParameterT
        );

    /**
     * @brief the index of the lower limit of the span in which the ParameterT lays.
     * @note From Piegl and Tiller, The NURBS Book, Algorithm A2.1
     */
    IndexType KRATOS_API(KRATOS_CORE) GetLowerSpan(
        const SizeType PolynomialDegree,
        const Vector& rKnots,
        const double ParameterT
        );

    /**
     * @brief Computes the degree of a nurbs/ b-spline shape by:
     * @param NumberOfKnots and
     * @param NumberOfControlPoints
     */
    SizeType KRATOS_API(KRATOS_CORE) GetPolynomialDegree(
        const SizeType NumberOfKnots, 
        const SizeType NumberOfControlPoints
        );
    
    /**
     * @brief Computes the number of knots of a nurbs/ b-spline shape by:
     * @param PolynomialDegree and
     * @param NumberOfControlPoints
     */
    SizeType KRATOS_API(KRATOS_CORE) GetNumberOfKnots(
        const SizeType PolynomialDegree, 
        const SizeType NumberOfControlPoints
        );

    /**
     * @brief Computes the number of control points of a nurbs/ b-spline shape by:
     * @param PolynomialDegree and
     * @param NumberOfKnots
     */
    SizeType KRATOS_API(KRATOS_CORE) GetNumberOfControlPoints(
        const SizeType PolynomialDegree, 
        const SizeType NumberOfKnots
        );

    /**
     * @brief Computes the number of spans of a nurbs/ b-spline shape by:
     * @param PolynomialDegree and
     * @param NumberOfKnots
     */
    SizeType KRATOS_API(KRATOS_CORE) GetNumberOfSpans(
        const SizeType PolynomialDegree, 
        const SizeType NumberOfKnots
        );

    /**
     * @brief Computes the binomial coefficient for (N || K).
     */
    constexpr inline SizeType GetBinomCoefficient(
        const SizeType N, 
        const SizeType K
        ) noexcept
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

    /**
     * @brief Computes a vector index from two matrix indicies.
     * @return index within vector
     */
    constexpr inline IndexType GetVectorIndexFromMatrixIndices(
        const SizeType NumberPerRow, 
        const SizeType NumberPerColumn,
        const IndexType RowIndex, 
        const IndexType ColumnIndex
        ) noexcept
    {
        return ColumnIndex * NumberPerRow + RowIndex;
    }


    /**
     * @brief Computes two matrix indices from vector index.
     * @return indices within Matrix
     */
    inline std::pair<IndexType, IndexType> KRATOS_API(KRATOS_CORE) GetMatrixIndicesFromVectorIndex(
        const SizeType NumberPerRow,
        const SizeType NumberPerColumn,
        const IndexType Index
        ) noexcept;
        
    ///@}
}; // class NurbsUtility
///@}
} // namespace Kratos

#endif // KRATOS_NURBS_UTILITY_H_INCLUDED defined 
