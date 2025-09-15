//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#pragma once

// System includes
#include <array>
#include <vector>

// External includes

// Project includes
#include "includes/kratos_export_api.h"

namespace Kratos::PolynomialUtilities
{
    using PolynomialType = std::vector<double>;
    using IntervalType = std::array<double, 2>;

    /** Degree of the polynomial */
    std::size_t KRATOS_API(KRATOS_CORE) Degree(const PolynomialType& rPolynomial);

    /** Evaluate the polynomial at the given x value */
    double KRATOS_API(KRATOS_CORE) Evaluate(const PolynomialType& rPolynomial, double x);

    /** Derivative of the polynomial */
    PolynomialType KRATOS_API(KRATOS_CORE) Differentiate(const PolynomialType& rPolynomial);

    /** Product of two polynomials */
    PolynomialType KRATOS_API(KRATOS_CORE) Multiply(const PolynomialType& rA, const PolynomialType& rB);

    /** Euclidian division of polynomials rA / rB = rQuotient*rB + rRemainder */
    void KRATOS_API(KRATOS_CORE) Divide(
        PolynomialType& rQuotient,
        PolynomialType& rRemainder,
        const PolynomialType& rA,
        const PolynomialType& rB
        );

    /** 
     * @brief Define disjoint subintervals of rRange, each containing a single root of rPolynomial.
     * @details The function uses the signs of the Sturm sequence in order to isolate the roots.
     * @param rRootIntervals each term defines a subinterval of rRange containing a single root.
     * @param rPolynomial the polynomial for which we want to find the roots.
     * @param rRange the interval in which the roots are to be found.
     */
    void KRATOS_API(KRATOS_CORE) IsolateRoots(
        std::vector<IntervalType>& rRootIntervals,
        const PolynomialType& rPolynomial,
        const IntervalType& rRange
        );

    /** 
     * @brief Find a root of rPolynomial within the interval defined by rRange.
     * @details The implementation assumes that rRange contains exactly one root
     * and that rRange[0]*rRange[1] < 0
     * @param rPolynomial the polynomial
     * @param rRange interval containing the root
     * @result The root.
     * @see IsolateRoots to identify suitable ranges.
     */
    double KRATOS_API(KRATOS_CORE) FindRoot(
        const PolynomialType& rPolynomial,
        const IntervalType& rRange
        );

}