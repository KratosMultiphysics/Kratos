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

#if !defined(KRATOS_NURBS_CURVE_TESSELLATION_H_INCLUDED )
#define  KRATOS_NURBS_CURVE_TESSELLATION_H_INCLUDED

#include "geometries/geometry.h"

#include "includes/ublas_interface.h"
#include "containers/array_1d.h"

namespace Kratos {

template <typename TCurve>
class NurbsCurveTessellation
{
public:     // types
    using Vector = array_1d<double, TCurve::GetDimension()>;
    using ParameterPoint = std::pair<double, Vector>;

private:    // static methods
    static double Cross(const array_1d<double, 2>& VectorU,
        const array_1d<double, 2>& VectorV)
    {
        return VectorV[0] * VectorU[1] - VectorV[1] * VectorU[0];
    }

    static array_1d<double, 3> Cross(const array_1d<double, 3>& VectorU,
        const array_1d<double, 3>& VectorV)
    {
        array_1d<double, 3> result;
        result[0] = VectorV[1] * VectorU[2] - VectorV[2] * VectorU[1];
        result[1] = VectorV[2] * VectorU[0] - VectorV[0] * VectorU[2];
        result[2] = VectorV[0] * VectorU[1] - VectorV[1] * VectorU[0];
        return result;
    }

    static double Norm(double Value)
    {
        return std::abs(Value);
    }

    static double Norm(const array_1d<double, 2>& Vector)
    {
        return std::sqrt(Vector[0] * Vector[0] + Vector[1] * Vector[1]);
    }

    static double Norm(const array_1d<double, 3>& Vector)
    {
        return std::sqrt(Vector[0] * Vector[0] + Vector[1] * Vector[1] +
            Vector[2] * Vector[2]);
    }

    static double DistanceToLine(const Vector& Point, const Vector& LineA,
        const Vector& LineB)
    {
        Vector vector_v = LineA - Point;
        Vector vector_u = LineB - LineA;

        return Norm(Cross(vector_v, vector_u)) / Norm(vector_u);
    }

public:     // static methods
    static std::vector<ParameterPoint> Compute(const TCurve& Curve,
        const Interval Domain, const double Tolerance)
    {
        std::vector<ParameterPoint> sample_points;
        std::vector<ParameterPoint> points;

        // compute sample points

        for (const auto& span : Curve.GetSpans()) {
            const Interval normalized_span = Domain.GetNormalizedInterval(span);

            if (normalized_span.GetLength() < 1e-7) {
                continue;
            }

            const double t = normalized_span.GetT0();
            const Vector point = Curve.GetPointAt(span.GetT0());

            sample_points.emplace_back(t, point);
        }

        sample_points.emplace_back(1.0,
            Curve.GetPointAt(Domain.GetParameterAtNormalized(1.0)));

        std::sort(std::begin(sample_points), std::end(sample_points),
            [](const ParameterPoint& lhs, const ParameterPoint& rhs) {
                return std::get<0>(lhs) > std::get<0>(rhs);
            }
        );

        // compute polyline

        const int n = Curve.GetDegree() * 2 + 1;

        while (true) {
            const auto parameter_point_a = sample_points.back();

            const auto t_a = std::get<0>(parameter_point_a);
            const auto point_a = std::get<1>(parameter_point_a);

            sample_points.pop_back();

            points.emplace_back(Domain.GetParameterAtNormalized(t_a), point_a);

            if (sample_points.size() == 0) {
                break;
            }

            while (true) {
                const auto parameter_point_b = sample_points.back();

                const auto t_b = std::get<0>(parameter_point_b);
                const auto point_b = std::get<1>(parameter_point_b);

                double max_distance {0};
                ParameterPoint max_point;

                for (int i = 1; i <= n; i++) {
                    const double t = Interval::GetParameterAtNormalized(t_a,
                        t_b, i / double(n + 1));
                    const Vector point = Curve.GetPointAt(
                        Domain.GetParameterAtNormalized(t));

                    const double distance = DistanceToLine(point, point_a,
                        point_b);

                    if (distance > max_distance) {
                        max_distance = distance;
                        max_point = {t, point};
                    }
                }

                if (max_distance < Tolerance) {
                    break;
                }

                sample_points.push_back(max_point);
            }
        }

        return points;
    }
}; 

} // namespace NurbsCurveTessellation

#endif // KRATOS_NURBS_CURVE_TESSELLATION_H_INCLUDED defined