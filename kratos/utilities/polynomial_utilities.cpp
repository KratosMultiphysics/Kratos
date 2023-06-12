//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela Dalmau
//

// System includes
#include <cmath>

// External includes

// Project includes
#include "polynomial_utilities.h"
#include "includes/define.h"
#include "utilities/brent_iteration.h"

namespace Kratos::PolynomialUtilities {

namespace {
    constexpr double DROP_TOLERANCE = 1e-12;
    constexpr double BRENT_TOLERANCE = 1e-12;
    constexpr int BRENT_MAX_ITER = 100;

    void DropLeadingZeros(PolynomialType& rPolynomial) {
        std::size_t offset = 0;
        std::size_t size = rPolynomial.size();

        // Note: size - 1 because the last term is never dropped,
        // the zero polynomial is {0.0}
        while (std::abs(rPolynomial[offset]) < DROP_TOLERANCE && offset < size-1) {
            ++offset;
        }

        if (offset == 0) return;

        for (std::size_t i=0; i < size - offset; i++) {
            rPolynomial[i] = rPolynomial[i+offset];
        }
        rPolynomial.resize(size-offset);
    }

    // Compute the next term in the Sturm series
    PolynomialType Sturm(const PolynomialType& rS1, const PolynomialType& rS2) {
        PolynomialType q, r;
        Divide(q, r, rS1, rS2);
        for (auto& c: r) {
            c *= -1;
        }
        return r;
    }

    std::size_t SignChanges(const std::vector<PolynomialType>& rSeries, double Coordinate) {
        double value, last = Evaluate(rSeries[0], Coordinate);
        std::size_t count = 0;

        for (auto iter = rSeries.begin() + 1; iter < rSeries.end(); ++iter) {
            value = Evaluate(*iter, Coordinate);
            if (last*value < 0) ++count;
            last = value;
        }

        return count;
    }
}

std::size_t Degree(const PolynomialType& rPolynomial)
{
    return rPolynomial.size() - 1;
}

double Evaluate(const PolynomialType& rPolynomial, double x)
{
    auto iter = rPolynomial.begin();
    double value = *(iter++);
    for (; iter != rPolynomial.end(); ++iter) {
        value = value*x + *iter;
    }

    return value;
}

PolynomialType Differentiate(const PolynomialType& rPolynomial)
{
    const auto degree = Degree(rPolynomial);
    PolynomialType deriv;
    deriv.reserve(rPolynomial.size() - 1);

    auto iter = rPolynomial.begin();
    for (std::size_t i = degree; i > 0; i--) {
        deriv.push_back(i*(*iter++));
    }

    return deriv;
}

PolynomialType Multiply(const PolynomialType& rA, const PolynomialType& rB)
{
    const auto size_a = rA.size();
    const auto size_b = rB.size();
    PolynomialType product(size_a + size_b - 1, 0.0);

    for (std::size_t i = 0; i < size_b; i++) {
        const double coeff = rB[i];
        for (std::size_t j = 0; j < size_a; j++) {
            product[i+j] += coeff*rA[j];
        }
    }

    return product;
}

void Divide(
    PolynomialType& rQuotient,
    PolynomialType& rRemainder,
    const PolynomialType& rA,
    const PolynomialType& rB
    )
{
    const std::size_t deg_a = Degree(rA);
    const std::size_t deg_b = Degree(rB);
    const std::size_t deg_q = deg_a - deg_b + 1;

    rQuotient.clear();
    rQuotient.reserve(deg_q);

    rRemainder = rA;

    // TODO: multiply actually multiplies by a few zeros here...
    // TODO: assuming that neither rA or rB have leading zeros
    for (std::size_t i = 0; i < deg_q; i++) {
        double s = rRemainder[0] / rB[0];
        rQuotient.push_back(s);
        PolynomialType aux(rRemainder.size() - deg_b, 0.0);
        aux[0] = s;
        aux = Multiply(aux, rB);
        // The leading term of remainder should be zero and is discarded
        // by writing all terms one position earlier and popping the last term
        for (std::size_t j = 1; j < rRemainder.size(); j++) {
            rRemainder[j-1] = rRemainder[j] - aux[j];
        }
        rRemainder.pop_back();
    }

    DropLeadingZeros(rRemainder);
}

void IsolateRoots(
    std::vector<IntervalType>& rRootIntervals,
    const PolynomialType& rPolynomial,
    const IntervalType& rRange
    )
{
    std::size_t sturm_seq_size = rPolynomial.size();
    std::vector<PolynomialType> sturm_sequence{rPolynomial, Differentiate(rPolynomial)};
    for (std::size_t i = 2; i < sturm_seq_size; i++) {
        sturm_sequence.push_back(Sturm(sturm_sequence[i-2], sturm_sequence[i-1]));
    }

    std::size_t va = SignChanges(sturm_sequence, rRange[0]);
    std::size_t vb = SignChanges(sturm_sequence, rRange[1]);
    std::size_t nroots = (vb > va) ? vb - va : va - vb;

    rRootIntervals.reserve(nroots);
    rRootIntervals.resize(0);

    if (nroots == 0) return;

    std::vector<IntervalType> candidates{rRange};
    using IntervalCountType = std::array<std::size_t,2>;
    std::vector<IntervalCountType> candidate_counts{{va, vb}};

    while (candidates.size() > 0) {
        IntervalType range = *(candidates.rbegin());
        auto counts = *(candidate_counts.rbegin());

        double c = 0.5*(range[0] + range[1]);
        std::size_t vc = SignChanges(sturm_sequence, c);

        candidates.pop_back();
        candidate_counts.pop_back();

        // roots in [a,c]
        auto a = range[0];
        auto va = counts[0];
        nroots = (vc > va) ? vc - va : va - vc;
        if (nroots == 1) {
            rRootIntervals.push_back(IntervalType{a, c});
        } else if (nroots > 1) {
            candidates.push_back(IntervalType{a, c});
            candidate_counts.push_back(IntervalCountType{va, vc});
        }

        // roots in [c,b]
        auto b = range[1];
        auto vb = counts[1];
        nroots = (vc > vb) ? vc - vb : vb - vc;
        if (nroots == 1) {
            rRootIntervals.push_back(IntervalType{c, b});
        } else if (nroots > 1) {
            candidates.push_back(IntervalType{c, b});
            candidate_counts.push_back(IntervalCountType{vc, vb});
        }
    }
}

double FindRoot(
    const PolynomialType& rPolynomial,
    const IntervalType& rRange
    )
{
    const auto func = [&rPolynomial](double Coordinate){
        return Evaluate(rPolynomial, Coordinate);
    };

    return BrentIteration::FindRoot(
        func, rRange[0], rRange[1], BRENT_TOLERANCE, BRENT_MAX_ITER);
}

}