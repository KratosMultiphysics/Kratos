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

#if !defined(KRATOS_PROJECTION_NURBS_GEOMETRY_UTILITIES_H_INCLUDED)
#define KRATOS_PROJECTION_NURBS_GEOMETRY_UTILITIES_H_INCLUDED

// System includes

// External includes

namespace Kratos
{
    namespace ProjectionNurbsGeometryUtilities
    {
        typedef array_1d<double, 3> CoordinatesArrayType;


        template <int TDimension>
        static bool NewtonRaphsonCurve(
            const double rInitialGuessParameter,
            const CoordinatesArrayType& rPoint
            CoordinatesArrayType& rResult,
            const NurbsCurveGeometry<TDimension>& rNurbsCurve,
            const int MaxIterations,
            const double ModelTolerance,
            const double Accuracy)
        {
            rResult[0] = rInitialGuessParameter;

            double min = rNurbsCurve.DomainInterval().Min();
            double max = rNurbsCurve.DomainInterval().Max();

            for (int i = 0; i < MaxIterations; i++) 
            {
                auto derivatives = rNurbsCurve->GlobalDerivatives(
                    rParameter,
                    2);

                array_1d<double, 3> distance_vector = derivatives[0] - rPoint;
                double distance = norm_2(distance_vector);

                double c2n = inner_prod(derivatives[1], distance_vector);
                double c2d = norm_2(derivatives[1]) * distance;
                double c2v = c2d != 0 
                    ? c2n / c2d
                    : 0;

                // Break condition
                if (distance < ModelTolerance
                    && std::abs(c2v) < Accuracy)
                    return true;

                double delta_t = inner_prod(derivatives[1], distance_vector)
                    / (inner_prod(derivatives[2], distance_vector) + pow(norm_2(derivatives[1]), 2));

                rParameter -= delta_t;

                //Alternative if (rNurbsCurve.DomainInterval().IsInside())
                if (rParameter < min || rParameter > max)
                {
                    return false;
                }
            }

            return false;
        }
    }

    template <int TDimension>
    static bool NewtonRaphsonSurface(
        const CoordinatesArrayType rInitialGuessParameter,
        const CoordinatesArrayType& rPoint
        CoordinatesArrayType& rResult,
        const NurbsSurfaceGeometry<TDimension>& rNurbsSurface,
        const int MaxIterations,
        const double ModelTolerance,
        const double Accuracy)
    {
        double rPoint[0] = rInitialGuessParameter[0];
        double rPoint[1] = rInitialGuessParameter[1];

        for (int i = 0; i < MaxIterations; i++) {
            const auto s = pSurface->GlobalDerivatives(rPoint[0], rPoint[1], 2);

            const array_1d<double, 3> distance = s[0] - rPoint;

            const double c1v = norm_2(distance);

            if (c1v < Accuracy) {
                rNewLocalCoordinates[0] = rPoint[0];
                rNewLocalCoordinates[1] = rPoint[1];

                return true;
            }

            double s1_l = norm_2(s[1]);
            double s2_l = norm_2(s[2]);

            double c2an = std::abs(inner_prod(s[1], distance));
            double c2ad = s1_l * c1v;

            double c2bn = std::abs(inner_prod(s[2], distance));
            double c2bd = s2_l * c1v;

            double c2av = c2an / c2ad;
            double c2bv = c2bn / c2bd;

            if (c2av < Tolerance && c2bv < Tolerance) {
                rNewLocalCoordinates[0] = rPoint[0];
                rNewLocalCoordinates[1] = rPoint[1];

                return true;
            }

            double f = inner_prod(s[1], distance);
            double g = inner_prod(s[2], distance);

            double J_00 = inner_prod(s[1], s[1]) + inner_prod(s[3], distance);
            double J_01 = inner_prod(s[1], s[2]) + inner_prod(s[4], distance);
            double J_11 = inner_prod(s[2], s[2]) + inner_prod(s[5], distance);

            double det_J = J_00 * J_11 - J_01 * J_01;

            double d_u = (g * J_01 - f * J_11) / det_J;
            double d_v = (f * J_01 - g * J_00) / det_J;

            rPoint[0] += d_u;
            rPoint[1] += d_v;
        }

        rNewLocalCoordinates[0] = rPoint[0];
        rNewLocalCoordinates[1] = rPoint[1];
        return false;
    }
}

} // namespace Kratos

#endif // KRATOS_PROJECTION_NURBS_GEOMETRY_UTILITIES_H_INCLUDED
