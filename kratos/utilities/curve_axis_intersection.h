//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//


#if !defined(KRATOS_CURVE_AXIS_INTERSECTION_INCLUDED )
#define  KRATOS_CURVE_AXIS_INTERSECTION_INCLUDED


// System includes

// External includes

// Project includes
#include "nurbs_curve_tessellation.h"

namespace Kratos
{
    template<std::size_t TWorkingSpaceDimension, class TNodeType>
    class CurveAxisIntersection
    {
    public:
        typedef std::size_t SizeType;
        typedef std::size_t IndexType;

        typedef Geometry<TNodeType> GeometryType;
        typedef typename GeometryType::Pointer GeometryPointerType;
        typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;

        typedef NurbsCurveTessellation<TWorkingSpaceDimension, PointerVector<TNodeType>> CurveTesselationType;

    private:
        static double BisectionToAxis(
            const GeometryType& rCurve,
            double IntersectionAxis1,
            double IntersectionAxis2,
            double Parameter1,
            double Parameter2,
            IndexType AxisDirectionIndex,
            double Tolerance = 1e-6)
        {
            // obtain correct axis
            double max = std::max(Parameter1, Parameter2);
            double min = std::min(Parameter1, Parameter2);
            const double intersection_axis = ((IntersectionAxis1 > min) && (IntersectionAxis1 < max))
                ? IntersectionAxis1
                : IntersectionAxis2;

            double parameter_smaller, parameter_bigger;
            CoordinatesArrayType point_1, point_2;
            CoordinatesArrayType parameter_1 = ZeroVector(3);
            parameter_1[0] = Parameter1;
            CoordinatesArrayType parameter_2 = ZeroVector(3);
            parameter_2[0] = Parameter2;
            rCurve.GlobalCoordinates(point_1, parameter_1);
            rCurve.GlobalCoordinates(point_2, parameter_2);

            double distance = point_1[AxisDirectionIndex] - intersection_axis;
            if (distance < 0) {
                parameter_smaller = parameter_1[0];
                parameter_bigger = parameter_2[0];
            }
            else {
                parameter_smaller = parameter_2[0];
                parameter_bigger = parameter_1[0];
            }

            for (IndexType i = 0; i <= IndexType(50); ++i)
            {
                CoordinatesArrayType new_parameter = ZeroVector(3);
                new_parameter[0] = (parameter_smaller + parameter_bigger) / 2;
                CoordinatesArrayType point_new;
                rCurve.GlobalCoordinates(point_new, new_parameter);
                distance = point_new[AxisDirectionIndex] - intersection_axis;
                if (std::abs(distance) < Tolerance)
                {
                    return new_parameter[0];
                }
                if (distance < 0)
                    parameter_smaller = new_parameter[0];
                else
                    parameter_bigger = new_parameter[0];
            }

            return 0.0;
        }

        static void GetSpanIndex(
            const std::vector<double>& rAxis,
            IndexType& rSpanIndex,
            double& rMin,
            double& rMax,
            double Parameter)
        {
            KRATOS_ERROR_IF(rAxis.size() < 2)
                << "Size of axis vector needs to be at least 2. It must contain the boundaries, too."
                << std::endl;

            for (IndexType i = 0; i < rAxis.size(); ++i) {
                KRATOS_ERROR_IF(i == rAxis.size() - 1)
                    << "Point of polygon not within the axis boundaries." << std::endl;
                rMin = std::min(rAxis[i], rAxis[i + 1]);
                rMax = std::max(rAxis[i], rAxis[i + 1]);

                if ((Parameter > rMin) && (Parameter < rMax)) {
                    rSpanIndex = i;
                    return;
                }
            }
        }

    public:
        static std::vector<double> ComputeAxisIntersection(
            const GeometryType& rGeometry,
            double Start,
            double End,
            const std::vector<double>& rAxis1,
            const std::vector<double>& rAxis2,
            const double Tolerance)
        {
            std::vector<double> intersection_parameters;

            intersection_parameters.push_back(Start);

            // approximate curve with a polyline

            const auto polygon = CurveTesselationType::ComputePolygon(
                rGeometry, 100, Start, End);

            // initialize axes
            IndexType axis_index_1, axis_index_2;
            double min_1, max_1, min_2, max_2;
            GetSpanIndex(rAxis1, axis_index_1, min_1, max_1, std::get<1>(polygon[0])[0]);
            GetSpanIndex(rAxis2, axis_index_2, min_2, max_2, std::get<1>(polygon[0])[1]);

            for (IndexType i = 1; i < polygon.size(); ++i) {
                if (std::get<1>(polygon[i])[0] < min_1 || std::get<1>(polygon[i])[0] > max_1) {
                    double intersection_parameter = BisectionToAxis(
                        rGeometry, rAxis1[axis_index_1], rAxis1[axis_index_1 + 1],
                        std::get<0>(polygon[i]), std::get<0>(polygon[i + 1]), Tolerance);
                    intersection_parameters.push_back(intersection_parameter);
                    GetSpanIndex(rAxis1, axis_index_1, min_1, max_1, std::get<1>(polygon[i])[0]);
                }

                if (std::get<1>(polygon[i])[1] < min_2 || std::get<1>(polygon[i])[1] > max_2) {
                    double intersection_parameter = BisectionToAxis(
                        rGeometry, rAxis2[axis_index_2], rAxis2[axis_index_2 + 1],
                        std::get<0>(polygon[i]), std::get<0>(polygon[i + 1]), Tolerance);
                    intersection_parameters.push_back(intersection_parameter);
                    GetSpanIndex(rAxis2, axis_index_2, min_2, max_2, std::get<1>(polygon[i])[1]);
                }
            }

            if (std::abs(intersection_parameters[intersection_parameters.size() - 1] - End) > Tolerance)
                intersection_parameters.push_back(End);

            std::sort(std::begin(intersection_parameters), std::end(intersection_parameters));

            auto last = std::unique(std::begin(intersection_parameters), std::end(intersection_parameters),
                [=](double a, double b) { return b - a < Tolerance; });

            auto nb_unique = std::distance(std::begin(intersection_parameters), last);

            intersection_parameters.resize(nb_unique);

            return intersection_parameters;
        }
    };
} // namespace Kratos.

#endif // KRATOS_CURVE_AXIS_INTERSECTION_INCLUDED  defined
