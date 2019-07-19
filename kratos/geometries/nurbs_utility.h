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

class NurbsUtility
{
public:     // static methods
    template <typename TKnots>
    static int GetUpperSpan(const int Degree, const TKnots& rKnots,
        const double ParameterT)
    {
        auto span = std::upper_bound(std::begin(rKnots) + Degree,
            std::end(rKnots) - Degree, ParameterT) - std::begin(rKnots) - 1;
        return static_cast<int>(span);
    }

    template <typename TKnots>
    static int GetLowerSpan(const int Degree, const TKnots& rKnots,
        const double ParameterT)
    {
        auto span = std::lower_bound(std::begin(rKnots) + Degree,
            std::end(rKnots) - Degree, ParameterT) - std::begin(rKnots) - 1;
        return static_cast<int>(span);
    }

    static int GetDegree(const int NbKnots, const int NbPoles)
    {
        return NbKnots - NbPoles + 1;
    }

    static int GetNbKnots(const int Degree, const int NbPoles)
    {
        return NbPoles + Degree - 1;
    }

    static int GetNbPoles(const int Degree, const int NbKnots)
    {
        return NbKnots - Degree + 1;
    }

    static int GetNbSpans(const int Degree, const int NbKnots)
    {
        return NbKnots - 2 * Degree + 1;
    }

    static constexpr inline int Binom(const int N, const int K) noexcept
    {
        // clang-format off
        return
            (K > N               ) ? 0           :  // out of range
            (K == 0 || K == N    ) ? 1           :  // edge
            (K == 1 || K == N - 1) ? N           :  // first
            (K + K < N           ) ?                // recursive:
            (Binom(N - 1, K - 1  ) * N) / K      :  //   path to K = 1     faster
            (Binom(N - 1, K      ) * N) / (N - K);  //   path to K = n - 1 faster
        // clang-format on
    }

    static constexpr inline int GetSingleIndex(const int Rows, const int Cols,
        const int Row, const int Col) noexcept
    {
        return Row * Cols + Col;
    }

    static constexpr inline std::pair<int, int> GetDoubleIndex(const int Rows,
        const int Cols, const int Index) noexcept
    {
        const int row = Index / Cols;
        const int col = Index % Cols;
        return {row, col};
    }
}; // class NurbsUtility

} // namespace Kratos

#endif // KRATOS_NURBS_UTILITY_H_INCLUDED defined