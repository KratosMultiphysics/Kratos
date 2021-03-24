//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Andreas Apostolatos
//					 Tobias Teschemacher
//					 Thomas Oberbichler
//					
//
//  Ported from the ANurbs library (https://github.com/oberbichler/ANurbs)
//

#if !defined(KRATOS_NURBS_CURVE_TESSELLATION_H_INCLUDED )
#define  KRATOS_NURBS_CURVE_TESSELLATION_H_INCLUDED

#include "utilities/math_utils.h"
#include "geometries/geometry.h"
#include "geometries/nurbs_curve_geometry.h"

#include "includes/ublas_interface.h"
#include "containers/array_1d.h"

namespace Kratos {

template <int TWorkingSpaceDimension, class TContainerPointType>
class NurbsCurveTessellation
{
public:

    ///@name Type Definitions
    ///@{

    typedef Geometry<typename TContainerPointType::value_type> GeometryType;
    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;
    typedef NurbsCurveGeometry<TWorkingSpaceDimension, TContainerPointType> NurbsCurveGeometryType;
    typedef std::vector<std::pair<double, Vector>> TessellationType;
    typedef typename GeometryType::IndexType IndexType;
    typedef typename GeometryType::SizeType SizeType;

    ///@}
private:
    ///@name Private Static Methods
    ///@{

    static double DistanceToLine(
        const typename GeometryType::CoordinatesArrayType& rPoint, 
        const typename GeometryType::CoordinatesArrayType& rLineA,
        const typename GeometryType::CoordinatesArrayType& rLineB
        )
    {
        typename GeometryType::CoordinatesArrayType vector_v = rLineA - rPoint;
        typename GeometryType::CoordinatesArrayType vector_u = rLineB - rLineA;

        return MathUtils<double>::Norm(MathUtils<double>::CrossProduct(vector_v, vector_u)) / MathUtils<double>::Norm(vector_u);
    }

    ///@}
public:
    ///@name Private Static Methods
    ///@{

    /// Conctructor for tessellation of a nurbs curve
    NurbsCurveTessellation()
    {
    }

    /** 
    * @brief This method tessellates a curve and stores the tessellation in the class
    * @param rGeometry Reference to the geometry
    * @param PolynomialDegree The polynomial degree of the curve
    * @param DomainInterval The curve interval which is to be tessellated
    * @param rKnotSpanIntervals Reference to the knot span intervals laying in the DomainInterval
    * @param Tolerance Tolerance for the choral error
    * @see ComputeTessellation
    */
    void Tessellate(
        const GeometryType& rGeometry, 
        const int PolynomialDegree,
        const NurbsInterval DomainInterval,
        const std::vector<NurbsInterval>& rKnotSpanIntervals,
        const double Tolerance)
    {
        mTesselation = ComputeTessellation<TWorkingSpaceDimension>(
            rGeometry,
            PolynomialDegree,
            DomainInterval,
            rKnotSpanIntervals,
            Tolerance);
    }

    /** 
    * @brief This method returns the tessellation of a curve
    * @param pGeometry Pointer to the geometry
    * @param PolynomialDegree The polynomial degree of the curve
    * @param DomainInterval The curve interval which is to be tessellated
    * @param KnotSpanIntervals The knot span intervals laying in the DomainInterval
    * @param Tolerance Tolerance for the choral error
    * @param ToSurfaceParameter defines if the tesselation is computed in
    *        global coordinates or in local coordinates of the underlying surface.
    * @return std::vector<std::pair<double, Vector>> tessellation
    * @see ANurbs library (https://github.com/oberbichler/ANurbs)
    */
    static TessellationType ComputeTessellation(
        const GeometryType& rGeometry,
        const int PolynomialDegree,
        const NurbsInterval DomainInterval,
        const std::vector<NurbsInterval>& rKnotSpanIntervals,
        const double Tolerance,
        bool ToSurfaceParameter = false
    )
    {
        TessellationType sample_points;
        TessellationType points;

        typename GeometryType::CoordinatesArrayType point;

        // compute sample points

        for (const auto& span : rKnotSpanIntervals) {
            const NurbsInterval normalized_span = DomainInterval.GetNormalizedInterval(span);

            if (normalized_span.GetLength() < 1e-7) {
                continue;
            }

            const double t = normalized_span.GetT0();
            typename GeometryType::CoordinatesArrayType t0;
            t0[0] = span.GetT0();

            ComputeGlobalCoordinates(
                point, t0, rGeometry, ToSurfaceParameter);

            sample_points.emplace_back(t, point);
        }

        typename GeometryType::CoordinatesArrayType t_at_normalized;
        t_at_normalized[0] = DomainInterval.GetParameterAtNormalized(1.0);

        ComputeGlobalCoordinates(
            point, t_at_normalized, rGeometry, ToSurfaceParameter);

        sample_points.emplace_back(1.0, point);

        std::sort(std::begin(sample_points), std::end(sample_points),
            [](std::pair<double, Vector> const& lhs, std::pair<double, Vector> const& rhs) {
                return std::get<0>(lhs) > std::get<0>(rhs);
            }
        );

        // compute polyline

        const int n = PolynomialDegree * 2 + 1;

        while (true) {
            const auto parameter_point_a = sample_points.back();

            const auto t_a = std::get<0>(parameter_point_a);
            const auto point_a = std::get<1>(parameter_point_a);

            sample_points.pop_back();

            points.emplace_back(DomainInterval.GetParameterAtNormalized(t_a), point_a);

            if (sample_points.size() == 0) {
                break;
            }

            while (true) {
                const auto parameter_point_b = sample_points.back();

                const auto t_b = std::get<0>(parameter_point_b);
                const auto point_b = std::get<1>(parameter_point_b);

                double max_distance{ 0 };
                std::pair<double, Vector> max_point;

                for (int i = 1; i <= n; i++) {
                    const double t = NurbsInterval::GetParameterAtNormalized(t_a,
                        t_b, i / (double(n) + 1.0));

                    t_at_normalized[0] = DomainInterval.GetParameterAtNormalized(t);

                    ComputeGlobalCoordinates(
                        point, t_at_normalized, rGeometry, ToSurfaceParameter);

                    const double distance = DistanceToLine(point, point_a,
                        point_b);

                    if (distance > max_distance) {
                        max_distance = distance;
                        max_point = { t, point };
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

    static void ComputeGlobalCoordinates(
        CoordinatesArrayType& rGlobalCoordinates,
        const CoordinatesArrayType& crLocaCoordinates,
        const GeometryType& rGeometry,
        bool to_surface_parameter = false
    )
    {
        if (!to_surface_parameter) {
            rGeometry.GlobalCoordinates(
                rGlobalCoordinates, crLocaCoordinates);
            return;
        }
        rGlobalCoordinates = crLocaCoordinates;
        rGeometry.Calculate(PARAMETER_2D_COORDINATES, rGlobalCoordinates);
    }

    /* @brief This method returns polygon of this curve with equal curve segments.
     * @param pGeometry Pointer to the geometry
     * @param NumberOfPoints The total amount of nodes including start and end node.
     * @param Start parameter of polygon.
     * @param End parameter of polygon.
     */
    static TessellationType ComputePolygon(
        const GeometryType& rGeometry,
        const SizeType NumberOfPoints,
        const double Start,
        const double End)
    {
        TessellationType points(NumberOfPoints);

        CoordinatesArrayType parameter = ZeroVector(3);

        double length = End - Start;
        double delta_length = length / (NumberOfPoints - 1);

        // compute sample points
        for (IndexType i = 0; i < NumberOfPoints; ++i) {
            parameter[0] = Start + delta_length * i;
            CoordinatesArrayType result;
            rGeometry.GlobalCoordinates(result, parameter);
            points[i] = { parameter[0], result };
        }

        return points;
    }

    /** 
    * @brief This method returns the already computed tessellation of a curve
    * @return return std::vector<std::pair<double, Vector>> tessellation
    */
    TessellationType GetTessellation() {
        return mTesselation;
    }

    private:
    ///@name Private Member Variables
    ///@{

    TessellationType mTesselation;

    ///@}
};

} // namespace NurbsCurveTessellation

#endif // KRATOS_NURBS_CURVE_TESSELLATION_H_INCLUDED defined
