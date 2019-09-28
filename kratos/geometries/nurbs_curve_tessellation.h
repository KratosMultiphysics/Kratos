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

// template <typename TCurve>
template <int TWorkingSpaceDimension, class TContainerPointType>
class NurbsCurveTessellation
{
public:

    ///@name Type Definitions
    ///@{

    typedef Geometry<typename TContainerPointType::value_type> GeometryType;
    typedef NurbsCurveGeometry<TWorkingSpaceDimension, TContainerPointType> NurbsCurveGeometryType;
    typedef typename GeometryType::IndexType IndexType;
    typedef typename GeometryType::SizeType SizeType;

private:    // static methods
    static double DistanceToLine(
        const typename GeometryType::CoordinatesArrayType& Point, 
        const typename GeometryType::CoordinatesArrayType& LineA,
        const typename GeometryType::CoordinatesArrayType& LineB
        )
    {
        typename GeometryType::CoordinatesArrayType vector_v = LineA - Point;
        typename GeometryType::CoordinatesArrayType vector_u = LineB - LineA;

        return MathUtils<double>::Norm(MathUtils<double>::CrossProduct(vector_v, vector_u)) / MathUtils<double>::Norm(vector_u);
    }

public:
    ///@}
    ///@name Life Cycle
    ///@{

    /// Conctructor for tessellation of a nurbs curve
    NurbsCurveTessellation()
    {
    }

    /** 
    * @brief This method tessellates a curve and stores the tessellation in the class
    * From ANurbs library (https://github.com/oberbichler/ANurbs)
    * @param pGeometry Pointer to the geometry
    * @param PolynomialDegree The polynomial degree of the curve
    * @param DomainInterval The curve interval which is to be tessellated
    * @param KnotSpanIntervals The knot span intervals laying in the DomainInterval
    * @param Tolerance Tolerance for the choral error
    * @see ComputeTessellation
    */
    void Tessellate(
        typename GeometryType::Pointer pGeometry, 
        int PolynomialDegree,
        Interval DomainInterval,
        std::vector<Interval> KnotSpanIntervals,
        const double Tolerance)
    {
        mTesselation = ComputeTessellation<TWorkingSpaceDimension>(
            pGeometry,
            PolynomialDegree,
            DomainInterval,
            KnotSpanIntervals,
            Tolerance);
    }

    /** 
    * @brief This method returns the tessellation of a curve
    * From ANurbs library (https://github.com/oberbichler/ANurbs)
    * @param pGeometry Pointer to the geometry
    * @param PolynomialDegree The polynomial degree of the curve
    * @param DomainInterval The curve interval which is to be tessellated
    * @param KnotSpanIntervals The knot span intervals laying in the DomainInterval
    * @param Tolerance Tolerance for the choral error
    * @return std::map<double, typename GeometryType::CoordinatesArrayType> with the coordinates in local and working space
    * @see ANurbs library (https://github.com/oberbichler/ANurbs)
    */
	static std::vector<std::pair<double, Vector>> ComputeTessellation(
        typename GeometryType::Pointer pGeometry,
        int PolynomialDegree,
        Interval DomainInterval,
        std::vector<Interval> KnotSpanIntervals,
        const double Tolerance
        )
    {
		static std::vector<std::pair<double, Vector>> sample_points;
		static std::vector<std::pair<double, Vector>> points;

        typename GeometryType::CoordinatesArrayType point;
        typename GeometryType::CoordinatesArrayType result;

        // compute sample points

        for (const auto& span : KnotSpanIntervals) {
            const Interval normalized_span = DomainInterval.GetNormalizedInterval(span);

            if (normalized_span.GetLength() < 1e-7) {
                continue;
            }

            const double t = normalized_span.GetT0();
            typename GeometryType::CoordinatesArrayType t0;
            t0[0] = span.GetT0();
            
            point = pGeometry->GlobalCoordinates(result, t0);

            sample_points.emplace_back(t, point);
        }

        typename GeometryType::CoordinatesArrayType t_at_normalized;
        t_at_normalized[0] = DomainInterval.GetParameterAtNormalized(1.0);

        point = pGeometry->GlobalCoordinates(result, t_at_normalized);

		sample_points.emplace_back(1.0, point);

		std::sort(std::begin(sample_points), std::end(sample_points),
			[](auto const& lhs, auto const& rhs) {
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

                double max_distance {0};
				std::pair<double, Vector> max_point;

                for (int i = 1; i <= n; i++) {
                    const double t = Interval::GetParameterAtNormalized(t_a,
                        t_b, i / (double(n) + 1.0));

                    t_at_normalized[0] = DomainInterval.GetParameterAtNormalized(t);

                    point = pGeometry->GlobalCoordinates(
                        result, t_at_normalized);

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

    private:
    ///@name Private Static Member Variables
    ///@{

    ///@}
    ///@name Private Member Variables
    ///@{

		std::vector<std::pair<double, Vector>> mTesselation;

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private Serialization
    ///@{
// 
    friend class Serializer;
};

} // namespace NurbsCurveTessellation

#endif // KRATOS_NURBS_CURVE_TESSELLATION_H_INCLUDED defined