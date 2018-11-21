#pragma once

#include "CurveGeometryBase.h"
#include "Pointer.h"

#include <stdexcept>
#include <vector>

namespace ANurbs {

template <typename TVector>
class CurveGeometry
    : public CurveGeometryBase<TVector>
{
public:
    using CurveGeometryBaseType = CurveGeometryBase<TVector>;
    using CurveGeometryType = CurveGeometry<TVector>;
    using typename CurveGeometryBaseType::KnotsType;
    using typename CurveGeometryBaseType::ScalarType;
    using typename CurveGeometryBaseType::VectorType;

protected:
    std::vector<VectorType> m_poles;
    std::vector<ScalarType> m_weights;

public:
    CurveGeometry(
        const int degree,
        const int nbPoles,
        const bool isRational)
        : CurveGeometryBaseType(degree, nbPoles)
        , m_poles(nbPoles)
        , m_weights(isRational ? nbPoles : 0)
    {
    }

    VectorType
    Pole(
        const int index) const override
    {
        return m_poles[index];
    }

    void
    SetPole(
        const int index,
        const VectorType& value) override
    {
        m_poles[index] = value;
    }

    bool
    IsRational() const override
    {
        return m_weights.size() != 0;
    }

    ScalarType
    Weight(
        const int index) const override
    {
        if (IsRational()) {
            return m_weights[index];
        } else {
            return 1;
        }
    }

    void
    SetWeight(
        const int index,
        const ScalarType value) override
    {
        if (IsRational()) {
            m_weights[index] = value;
        } else {
            throw std::invalid_argument("Geometry is not rational");
        }
    }

    Unique<CurveGeometryType>
    Refined(
        const std::vector<ScalarType>& knotsToInsert) const
    {
        int nbKnotsToInsert = static_cast<int>(knotsToInsert.size());

        int a = Knots::UpperSpan(this->Degree(), this->Knots(),
            knotsToInsert.front());
        int b = Knots::UpperSpan(this->Degree(), this->Knots(),
            knotsToInsert.back());

        int nbPolesRefined = this->NbPoles() + nbKnotsToInsert;
        int nbKnotsRefined = this->NbKnots() + nbKnotsToInsert;

        auto refined = New<CurveGeometryType>(this->Degree(), nbPolesRefined,
            true);

        for (int j = 0; j < a - this->Degree() + 2; j++) {
            refined->SetPole(j, this->WeightedPole(j));
            refined->SetWeight(j, Weight(j));
        }

        for (int j = b; j < this->NbPoles(); j++) {
            refined->SetPole(nbKnotsToInsert + j, this->WeightedPole(j));
            refined->SetWeight(nbKnotsToInsert + j, Weight(j));
        }

        for (int j = 0; j < a + 1; j++) {
            refined->SetKnot(j, this->Knot(j));
        }

        for (int j = b + this->Degree(); j < this->NbKnots(); j++) {
            refined->SetKnot(nbKnotsToInsert + j, this->Knot(j));
        }

        int i = b + this->Degree() - 1;
        int k = b + this->Degree() + nbKnotsToInsert - 1;

        for (int j = nbKnotsToInsert - 1; j > -1; j--) {
            while (knotsToInsert[j] <= this->Knot(i) && i > a) {
                refined->SetPole(k - this->Degree(),
                    this->WeightedPole(i - this->Degree()));
                refined->SetWeight(k - this->Degree(), Weight(i - this->Degree()));
                refined->SetKnot(k, this->Knot(i));
                k--;
                i--;
            }

            refined->SetPole(k - this->Degree(), refined->Pole(k -
                this->Degree() + 1));
            refined->SetWeight(k - this->Degree(), refined->Weight(k -
                this->Degree() + 1));

            for (int l = 1; l <= this->Degree(); l++) {
                int index = k - this->Degree() + l;
                ScalarType alpha = refined->Knot(k + l) - knotsToInsert[j];

                if (std::abs(alpha) == 0) {
                    refined->SetPole(index, refined->Pole(index + 1));
                    refined->SetWeight(index, refined->Weight(index + 1));
                } else {
                    alpha = alpha / (refined->Knot(k + l) - this->Knot(i + l -
                        this->Degree()));
                    refined->SetPole(index, refined->Pole(index) * alpha +
                        refined->Pole(index + 1) * (1 - alpha));
                    refined->SetWeight(index, refined->Weight(index) * alpha +
                        refined->Weight(index + 1) * (1 - alpha));
                }
            }

            refined->SetKnot(k, knotsToInsert[j]);

            k--;
        }

        for (int i = 0; i < nbPolesRefined; i++) {
            refined->SetPole(i, refined->Pole(i) * (1 / refined->Weight(i)));
        }

        return refined;
    }
};

using CurveGeometry1D = CurveGeometry<Point1D>;
using CurveGeometry2D = CurveGeometry<Point2D>;
using CurveGeometry3D = CurveGeometry<Point3D>;

} // namespace ANurbs
