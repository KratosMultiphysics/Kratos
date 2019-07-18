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
//                   Tobias Teschemacher
//

#if !defined(KRATOS_PROJECTION_ON_NURBS_CURVE_UTILITIES_H_INCLUDED)
#define KRATOS_PROJECTION_ON_NURBS_CURVE_UTILITIES_H_INCLUDED

// System includes

// External includes

namespace Kratos
{

    namespace ProjectionOnNurbsCurveUtilities
    {
        static bool ClosestPointToLine(
            double& rPointParameter,
            const array_1d<double, 3>& rPoint,
            const array_1d<double, 3>& rStartPoint,
            const array_1d<double, 3>& rEndPoint,
            const double& rStartParameter,
            const double& rEndParameter,
            const double& rAccuracy
        )
        {
            array_1d<double, 3> segment_vector = rEndPoint - rStartPoint;
            double segment_length = norm_2(segment_vector);

            array_1d<double, 3> segment_normal = segment_vector * (1.0 / segment_length);
            array_1d<double, 3> point_start_vector = rPoint - rStartPoint;
            double dot_product = inner_prod(point_start_vector, segment_normal);

            if (dot_product < 0) {
                rPointParameter = rStartParameter;
                return false;
            }

            if (dot_product > segment_length) {
                rPointParameter = rEndParameter;
                return false;
            }

            rPointParameter = rStartParameter + (rEndParameter - rStartParameter) * dot_product / segment_length;
            return true;
        }

        //template <int TDimension>
        static bool NewtonRaphson(
            double& rParameter,
            const double rInitialGuessParameter,
            array_1d<double, 3>& rPoint,
            const std::shared_ptr<CurveOnSurface<3>>& pCurve,
            const int MaxIterations,
            const double ModelTolerance,
            const double Accuracy)
        {
            rParameter = rInitialGuessParameter;

            double min = pCurve->Domain().Min();
            double max = pCurve->Domain().Max();

            for (int i = 0; i < MaxIterations; i++) 
            {
                auto derivatives = pCurve->DerivativesAt(rParameter, 2);

                array_1d<double, 3> distance_vector = derivatives[0] - rPoint;
                double distance = norm_2(distance_vector);

                double c2n = inner_prod(derivatives[1], distance_vector);
                double c2d = norm_2(derivatives[1]) * distance;
                double c2v = c2d != 0 
                    ? c2n / c2d
                    : 0;

                if (distance < ModelTolerance
                    && std::abs(c2v) < Accuracy)
                    return true;

                double delta_t = inner_prod(derivatives[1], distance_vector)
                    / (inner_prod(derivatives[2], distance_vector) + pow(norm_2(derivatives[1]), 2));

                rParameter -= delta_t;

                if (rParameter < min || rParameter > max)
                {
                    return false;
                }
            }

            return false;
        }
    }

} // namespace Kratos

#endif // KRATOS_PROJECTION_ON_NURBS_CURVE_UTILITIES_H_INCLUDED
