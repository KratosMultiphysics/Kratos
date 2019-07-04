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
    static int upper_span(const int degree, const TKnots& knots,
        const double& t)
    {
        auto span = std::upper_bound(std::begin(knots) + degree, std::end(knots)
            - degree, t) - std::begin(knots) - 1;
        return static_cast<int>(span);
    }

    template <typename TKnots>
    static int lower_span(const int degree, const TKnots& knots,
        const double& t)
    {
        auto span = std::lower_bound(std::begin(knots) + degree, std::end(knots)
            - degree, t) - std::begin(knots) - 1;
        return static_cast<int>(span);
    }

    static int degree(const int nb_knots, const int nb_poles)
    {
        return nb_knots - nb_poles + 1;
    }

    static int nb_knots(const int degree, const int nb_poles)
    {
        return nb_poles + degree - 1;
    }

    static int nb_poles(const int degree, const int nb_knots)
    {
        return nb_knots - degree + 1;
    }

    static int nb_spans(const int degree, const int nb_knots)
    {
        return nb_knots - 2 * degree + 1;
    }

    static constexpr inline int binom(const int n, const int k) noexcept
    {
        // clang-format off
        return
            (k > n               ) ? 0           :  // out of range
            (k == 0 || k == n    ) ? 1           :  // edge
            (k == 1 || k == n - 1) ? n           :  // first
            (k + k < n           ) ?                // recursive:
            (binom(n - 1, k - 1  ) * n) / k      :  //   path to k = 1     faster
            (binom(n - 1, k      ) * n) / (n - k);  //   path to k = n - 1 faster
        // clang-format on
    }

    static constexpr inline int single_index(const int rows, const int cols,
        const int row, const int col) noexcept
    {
        return row * cols + col;
    }

    static constexpr inline std::pair<int, int> double_index(const int rows,
        const int cols, const int index) noexcept
    {
        const int row = index / cols;
        const int col = index % cols;
        return {row, col};
    }
}; // class NurbsUtility

} // namespace Kratos

#endif // KRATOS_NURBS_UTILITY_H_INCLUDED defined