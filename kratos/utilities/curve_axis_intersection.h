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
#include <limits>

// External includes

// Project includes
#include "tessellation_utilities/curve_tessellation.h"

namespace Kratos
{
    template<class TNodeType>
    class CurveAxisIntersection
    {
    public:
        typedef std::size_t SizeType;
        typedef std::size_t IndexType;

        typedef Geometry<TNodeType> GeometryType;
        typedef typename GeometryType::Pointer GeometryPointerType;
        typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;

        typedef CurveTessellation<PointerVector<TNodeType>> CurveTesselationType;

    private:
        static double BisectionToAxis(
            const GeometryType& rCurve,
            double IntersectionAxis,
            double Parameter1,
            double Parameter2,
            IndexType AxisDirectionIndex,
            double Tolerance = 1e-6)
        {
            double parameter_smaller, parameter_bigger;
            CoordinatesArrayType point_1, point_2;
            CoordinatesArrayType parameter_1 = ZeroVector(3);
            parameter_1[0] = Parameter1;
            CoordinatesArrayType parameter_2 = ZeroVector(3);
            parameter_2[0] = Parameter2;
            rCurve.GlobalCoordinates(point_1, parameter_1);
            rCurve.GlobalCoordinates(point_2, parameter_2);

            double distance = point_1[AxisDirectionIndex] - IntersectionAxis;
            if (distance < -1e-8) {
                parameter_smaller = parameter_1[0];
                parameter_bigger = parameter_2[0];
            }
            else if (distance > 1e-8) {
                parameter_smaller = parameter_2[0];
                parameter_bigger = parameter_1[0];
            }
            else {
                if (norm_2(point_1) - norm_2(point_2) < 0)
                {
                    parameter_smaller = parameter_1[0];
                    parameter_bigger = parameter_2[0];
                }
                else{
                    parameter_smaller = parameter_2[0];
                    parameter_bigger = parameter_1[0];
                }
            }

            for (IndexType i = 0; i <= IndexType(50); ++i)
            {
                CoordinatesArrayType new_parameter = ZeroVector(3);
                new_parameter[0] = (parameter_smaller + parameter_bigger) / 2;
                CoordinatesArrayType point_new;
                rCurve.GlobalCoordinates(point_new, new_parameter);
                distance = point_new[AxisDirectionIndex] - IntersectionAxis;

                if (std::abs(distance) < Tolerance) {
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
            double Parameter,
            bool Ascending = true)
        {
            if (Ascending) {
                for (IndexType i = 0; i < rAxis.size() - 1; ++i) {
                    KRATOS_DEBUG_ERROR_IF(i == rAxis.size() - 1)
                        << "Point of polygon not within the axis boundaries. Axis are: "
                        << rAxis << ". Searched parameter is: " << Parameter << std::endl;

                    rMin = std::min(rAxis[i], rAxis[i + 1]);
                    rMax = std::max(rAxis[i], rAxis[i + 1]);

                    if ((Parameter >= rMin) && (Parameter < rMax)) {
                        rSpanIndex = i;
                        return;
                    }
                }
            }
            else {
                for (IndexType i = 0; i < rAxis.size() - 1; ++i) {
                    KRATOS_DEBUG_ERROR_IF(i == rAxis.size() - 1)
                        << "Point of polygon not within the axis boundaries. Axis are: "
                        << rAxis << ". Searched parameter is: " << Parameter << std::endl;

                    rMin = std::max(rAxis[i], rAxis[i + 1]);
                    rMax = std::min(rAxis[i], rAxis[i + 1]);

                    if ((Parameter >= rMin) && (Parameter < rMax)) {
                        rSpanIndex = i;
                        return;
                    }
                }
            }
        }

        /* @brief sort a std::vector<double> and delete duplicated entries.
         * @param resulting list to be sorted.
         * @param Tolerance between duplicated entries.
         */
        static void SortUnique(
            std::vector<double>& rIntersectionParameters,
            const double Tolerance)
        {
            std::sort(std::begin(rIntersectionParameters), std::end(rIntersectionParameters));

            auto last = std::unique(std::begin(rIntersectionParameters), std::end(rIntersectionParameters),
                [=](double a, double b) { return b - a < Tolerance; });

            auto nb_unique = std::distance(std::begin(rIntersectionParameters), last);

            rIntersectionParameters.resize(nb_unique);
        }

    public:
        static void ComputeAxisIntersection(
            std::vector<double>& rIntersectionParameters,
            const GeometryType& rGeometry,
            double Start,
            double End,
            const std::vector<double>& rAxis1,
            const std::vector<double>& rAxis2,
            const double Tolerance)
        {
            rIntersectionParameters.clear();
            rIntersectionParameters.push_back(Start);

            // linearise polygon
            const auto polygon = CurveTesselationType::ComputePolygon(
                rGeometry, 100, Start, End); //500

            KRATOS_ERROR_IF(rAxis1.size() < 2)
                << "Size of axis vector 1 is: " << rAxis1.size() << ". It needs to be at least of size 2. "
                << "It must contain the boundaries, too. Axis vector: " << rAxis1
                << std::endl;
            KRATOS_ERROR_IF(rAxis2.size() < 2)
                << "Size of axis vector 2 is: " << rAxis2.size() << ". It needs to be at least of size 2. "
                << "It must contain the boundaries, too. Axis vector: " << rAxis1
                << std::endl;

            bool ascending_1 = (rAxis1[0] < rAxis1[1]);
            bool ascending_2 = (rAxis2[0] < rAxis2[1]);

            // initialize axes
            IndexType axis_index_1, axis_index_2;
            double min_1 = std::numeric_limits<double>::max();
            double max_1 = std::numeric_limits<double>::lowest();
            double min_2 = std::numeric_limits<double>::max();
            double max_2 = std::numeric_limits<double>::lowest();
            GetSpanIndex(rAxis1, axis_index_1, min_1, max_1, std::get<1>(polygon[0])[0], ascending_1);
            GetSpanIndex(rAxis2, axis_index_2, min_2, max_2, std::get<1>(polygon[0])[1], ascending_2);
            double toll = 1e-13;

            // iterate through polygon and check for knot intersections
            for (IndexType i = 1; i < polygon.size(); ++i) {
                if (std::get<1>(polygon[i])[0] < min_1 - Tolerance) {
                    double intersection_parameter = BisectionToAxis(
                        rGeometry, min_1,
                        std::get<0>(polygon[i - 1]), std::get<0>(polygon[i]), 0, Tolerance);
                    rIntersectionParameters.push_back(intersection_parameter+toll);
                    GetSpanIndex(rAxis1, axis_index_1, min_1, max_1, std::get<1>(polygon[i])[0], ascending_1);
                }
                else if (std::get<1>(polygon[i])[0] > max_1 + Tolerance) {
                    double intersection_parameter = BisectionToAxis(
                        rGeometry, max_1,
                        std::get<0>(polygon[i - 1]), std::get<0>(polygon[i]), 0, Tolerance);
                    rIntersectionParameters.push_back(intersection_parameter+toll);
                    GetSpanIndex(rAxis1, axis_index_1, min_1, max_1, std::get<1>(polygon[i])[0], ascending_1);
                }

                if (std::get<1>(polygon[i])[1] < min_2 - Tolerance) {
                    double intersection_parameter = BisectionToAxis(
                        rGeometry, min_2,
                        std::get<0>(polygon[i - 1]), std::get<0>(polygon[i]), 1, Tolerance);
                    rIntersectionParameters.push_back(intersection_parameter+toll);
                    GetSpanIndex(rAxis2, axis_index_2, min_2, max_2, std::get<1>(polygon[i])[1], ascending_2);
                }
                else if (std::get<1>(polygon[i])[1] > max_2 + Tolerance) {
                    double intersection_parameter = BisectionToAxis(
                        rGeometry, max_2,
                        std::get<0>(polygon[i - 1]), std::get<0>(polygon[i]), 1, Tolerance);
                    rIntersectionParameters.push_back(intersection_parameter+toll);
                    GetSpanIndex(rAxis2, axis_index_2, min_2, max_2, std::get<1>(polygon[i])[1], ascending_2);
                }
            }

            // Add last element if different from knot intersection.
            if (std::abs(rIntersectionParameters[rIntersectionParameters.size() - 1] - End) > Tolerance)
                rIntersectionParameters.push_back(End);

            // sort and delete duplicated entries
            SortUnique(rIntersectionParameters, Tolerance);
        }
    };
} // namespace Kratos.

#endif // KRATOS_CURVE_AXIS_INTERSECTION_INCLUDED  defined
