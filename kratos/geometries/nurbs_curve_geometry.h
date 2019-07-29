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

#if !defined(KRATOS_NURBS_CURVE_GEOMETRY_H_INCLUDED )
#define  KRATOS_NURBS_CURVE_GEOMETRY_H_INCLUDED

#include "nurbs_curve_shape_functions.h"
#include "nurbs_interval.h"
#include "includes/ublas_interface.h"
#include "containers/array_1d.h"

#include <stdexcept>
#include <vector>

namespace Kratos {

template <int TDimension, typename TPoleContainer =
    std::vector<array_1d<double, TDimension>>>
class NurbsCurveGeometry
{
public:     // types
    using Type = NurbsCurveGeometry<TDimension>;
    using Vector = array_1d<double, TDimension>;

private:    // variables
    const int mDegree;
    std::vector<double> mKnots;
    TPoleContainer mPoles;
    std::vector<double> mWeights;

public:     // constructors
    NurbsCurveGeometry(
        const int Degree,
        int NbPoles,
        bool IsRational
    ) : mDegree(Degree),
        mPoles(NbPoles),
        mWeights(IsRational ? NbPoles : 0),
        mKnots(NurbsUtility::GetNbKnots(Degree, NbPoles))
    {
        static_assert(TDimension > 0, "Invalid dimension");
    }

    NurbsCurveGeometry(
        const int Degree,
        const std::vector<double>& rKnots,
        const std::vector<Vector>& rPoles
    ) : mDegree(Degree),
        mKnots(rKnots),
        mPoles(rPoles),
        mWeights()
    {
        static_assert(TDimension > 0, "Invalid dimension");

        if (rKnots.size() != NurbsUtility::GetNbKnots(Degree, rPoles.size())) {
            throw std::runtime_error("Number of knots and poles do not match");
        }
    }

    NurbsCurveGeometry(
        const int Degree,
        const std::vector<double>& rKnots,
        const std::vector<Vector>& rPoles,
        const std::vector<double>& rWeights
    ) : mDegree(Degree),
        mKnots(rKnots),
        mPoles(rPoles),
        mWeights(rWeights)
    {
        static_assert(TDimension > 0);

        if (rKnots.size() != NurbsUtility::GetNbKnots(Degree, rPoles.size())) {
            throw std::runtime_error("Number of knots and poles do not match");
        }

        if (rWeights.size() != rPoles.size()) {
            throw std::runtime_error(
                "Number of poles and weights do not match");
        }
    }

public:     // static methods
    static constexpr int GetDimension()
    {
        return TDimension;
    }

public:     // methods
    int GetDegree() const
    {
        return mDegree;
    }

    std::vector<Vector> GetDerivativesAt(const double ParameterT,
        const int Order) const
    {
        NurbsCurveShapeFunction shape_function;

        shape_function.Resize(mDegree, Order);

        if (mWeights.size() > 0) {
            shape_function.Compute(mKnots, [&](int i) { return GetWeight(i); }, ParameterT);
        } else {
            shape_function.Compute(mKnots, ParameterT);
        }

        std::vector<Vector> derivatives(shape_function.GetNbShapeFunctions());

        for (int order = 0; order < shape_function.GetNbShapeFunctions(); order++) {
            for (int i = 0; i < shape_function.GetNbNonzeroPoles(); i++) {
                int index = shape_function.GetFirstNonzeroPole() + i;

                if (i == 0) {
                    derivatives[order] = GetPole(index) * shape_function(order, i);
                } else {
                    derivatives[order] += GetPole(index) * shape_function(order, i);
                }
            }
        }

        return derivatives;
    }

    Interval GetDomain() const
    {
        return Interval(mKnots[mDegree - 1], mKnots[GetNbKnots() - mDegree]);
    }

    bool IsRational() const
    {
        return mWeights.size() != 0;
    }

    double GetKnot(const int Index) const
    {
        return mKnots[Index];
    }

    const std::vector<double>& GetKnots() const
    {
        return mKnots;
    }

    int GetNbKnots() const
    {
        return static_cast<int>(mKnots.size());
    }

    int GetNbPoles() const
    {
        return static_cast<int>(mPoles.size());
    }

    Vector GetPointAt(const double ParameterT) const
    {
        NurbsCurveShapeFunction shape_function;

        shape_function.Resize(mDegree, 0);

        if (mWeights.size() > 0) {
            shape_function.Compute(mKnots,
                [&](int i) { return GetWeight(i); }, ParameterT);
        } else {
            shape_function.Compute(mKnots, ParameterT);
        }

        Vector point;

        for (int i = 0; i < shape_function.GetNbNonzeroPoles(); i++) {
            const int index = shape_function.GetFirstNonzeroPole() + i;

            if (i == 0) {
                point = GetPole(index) * shape_function(0, i);
            } else {
                point += GetPole(index) * shape_function(0, i);
            }
        }

        return point;
    }

    Vector GetPole(const int Index) const
    {
        return mPoles.at(Index);
    }

    const TPoleContainer& GetPoles() const
    {
        return mPoles;
    }
    
    void SetKnot(const int Index, const double Value)
    {
        mKnots.at(Index) = Value;
    }
    
    void SetPole(const int Index, const Vector& rValue)
    {
        mPoles.at(Index) = rValue;
    }
    
    void SetWeight(const int Index, const double Value)
    {
        mWeights.at(Index) = Value;
    }

    std::pair<std::vector<int>, Matrix> GetShapeFunctionsAt(
        const double ParameterT, const int Order) const
    {
        NurbsCurveShapeFunction shape_function(GetDegree(), Order);

        if (IsRational()) {
            shape_function.Compute(GetKnots(), [&](int i) {
                return GetWeight(i); }, ParameterT);
        } else {
            shape_function.Compute(GetKnots(), ParameterT);
        }

        Matrix values(shape_function.GetNbShapeFunctions(),
            shape_function.GetNbNonzeroPoles());

        for (int i = 0; i < shape_function.GetNbShapeFunctions(); i++) {
            for (int j = 0; j < shape_function.GetNbNonzeroPoles(); j++) {
                values(i, j) = shape_function(i, j);
            }
        }

        return {shape_function.GetNonzeroPoleIndices(), values};
    }

    std::vector<Interval> Spans() const
    {
        const int first_span = GetDegree() - 1;
        const int last_span = GetNbKnots() - GetDegree() - 1;

        const int nb_spans = last_span - first_span + 1;

        std::vector<Interval> result(nb_spans);

        for (int i = 0; i < nb_spans; i++) {
            const double t0 = GetKnot(first_span + i);
            const double t1 = GetKnot(first_span + i + 1);

            result[i] = Interval(t0, t1);
        }

        return result;
    }
    
    double GetWeight(const int Index) const
    {
        return mWeights.at(Index);
    }
    
    const std::vector<double>& GetWeights() const
    {
        return mWeights;
    }
}; // class NurbsCurveGeometry

} // namespace Kratos

#endif // KRATOS_NURBS_CURVE_GEOMETRY_H_INCLUDED defined