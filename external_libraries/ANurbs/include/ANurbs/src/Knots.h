#pragma once

#include <algorithm>

namespace ANurbs {

class Knots
{
public:
    template <typename TScalar, typename TKnots>
    static int
    UpperSpan(
        const int degree,
        const TKnots& knots,
        const TScalar& t)
    {
        size_t span = std::upper_bound(std::begin(knots) + degree,
            std::end(knots) - degree, t) - std::begin(knots) - 1;
        return static_cast<int>(span);
    }

    template <typename TScalar, typename TKnots>
    static int
    LowerSpan(
        const int degree,
        const TKnots& knots,
        const TScalar& t)
    {
        size_t span = std::lower_bound(std::begin(knots) + degree,
            std::end(knots) - degree, t) - std::begin(knots) - 1;
        return static_cast<int>(span);
    }

    static int
    Degree(
        const int nbKnots,
        const int nbPoles)
    {
        return nbKnots - nbPoles + 1;
    }

    static int
    NbKnots(
        const int degree,
        const int nbPoles)
    {
        return nbPoles + degree - 1;
    }

    static int
    NbPoles(
        const int degree,
        const int nbKnots)
    {
        return nbKnots - degree + 1;
    }

    static int
    NbSpans(
        const int degree,
        const int nbKnots)
    {
        return nbKnots - 2 * degree + 1;
    }
};

} // namespace ANurbs
