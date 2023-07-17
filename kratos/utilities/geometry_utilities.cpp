//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

// System includes

// External includes

// Project includes
#include "utilities/geometry_utilities.h"
#include "utilities/geometrical_projection_utilities.h"

namespace Kratos
{

std::string GeometryUtils::GetGeometryName(const GeometryData::KratosGeometryType TypeOfGeometry)
{
    KRATOS_TRY;

    // Using switch over map as the compiler warns if some enum values are not handled in the switch
    switch(TypeOfGeometry) {
        case GeometryData::KratosGeometryType::Kratos_generic_type:
            return "GenericType";
        case GeometryData::KratosGeometryType::Kratos_Hexahedra3D20:
            return "Hexahedra3D20";
        case GeometryData::KratosGeometryType::Kratos_Hexahedra3D27:
            return "Hexahedra3D27";
        case GeometryData::KratosGeometryType::Kratos_Hexahedra3D8:
            return "Hexahedra3D8";
        case GeometryData::KratosGeometryType::Kratos_Prism3D15:
            return "Prism3D15";
        case GeometryData::KratosGeometryType::Kratos_Prism3D6:
            return "Prism3D6";
        case GeometryData::KratosGeometryType::Kratos_Pyramid3D13:
            return "Pyramid3D13";
        case GeometryData::KratosGeometryType::Kratos_Pyramid3D5:
            return "Pyramid3D5";
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4:
            return "Quadrilateral2D4";
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8:
            return "Quadrilateral2D8";
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9:
            return "Quadrilateral2D9";
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4:
            return "Quadrilateral3D4";
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8:
            return "Quadrilateral3D8";
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9:
            return "Quadrilateral3D9";
        case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10:
            return "Tetrahedra3D10";
        case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
            return "Tetrahedra3D4";
        case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
            return "Triangle2D3";
        case GeometryData::KratosGeometryType::Kratos_Triangle2D6:
            return "Triangle2D6";
        case GeometryData::KratosGeometryType::Kratos_Triangle2D10:
            return "Triangle2D10";
        case GeometryData::KratosGeometryType::Kratos_Triangle2D15:
            return "Triangle2D15";
        case GeometryData::KratosGeometryType::Kratos_Triangle3D3:
            return "Triangle3D3";
        case GeometryData::KratosGeometryType::Kratos_Triangle3D6:
            return "Triangle3D6";
        case GeometryData::KratosGeometryType::Kratos_Line2D2:
            return "Line2D2";
        case GeometryData::KratosGeometryType::Kratos_Line2D3:
            return "Line2D3";
        case GeometryData::KratosGeometryType::Kratos_Line2D4:
            return "Line2D4";
        case GeometryData::KratosGeometryType::Kratos_Line2D5:
            return "Line2D5";
        case GeometryData::KratosGeometryType::Kratos_Line3D2:
            return "Line3D2";
        case GeometryData::KratosGeometryType::Kratos_Line3D3:
            return "Line3D3";
        case GeometryData::KratosGeometryType::Kratos_Point2D:
            return "Point2D";
        case GeometryData::KratosGeometryType::Kratos_Point3D:
            return "Point3D";
        case GeometryData::KratosGeometryType::Kratos_Sphere3D1:
            return "Sphere3D1";
        case GeometryData::KratosGeometryType::Kratos_Nurbs_Curve:
            return "Nurbs_Curve";
        case GeometryData::KratosGeometryType::Kratos_Nurbs_Surface:
            return "Nurbs_Surface";
        case GeometryData::KratosGeometryType::Kratos_Nurbs_Volume:
            return "Nurbs_Volume";
        case GeometryData::KratosGeometryType::Kratos_Nurbs_Curve_On_Surface:
            return "Nurbs_Curve_On_Surface";
        case GeometryData::KratosGeometryType::Kratos_Surface_In_Nurbs_Volume:
            return "Surface_In_Nurbs_Volume";
        case GeometryData::KratosGeometryType::Kratos_Brep_Curve:
            return "Brep_Curve";
        case GeometryData::KratosGeometryType::Kratos_Brep_Surface:
            return "Brep_Surface";
        case GeometryData::KratosGeometryType::Kratos_Brep_Curve_On_Surface:
            return "Brep_Curve_On_Surface";
        case GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Geometry:
            return "Quadrature_Point_Geometry";
        case GeometryData::KratosGeometryType::Kratos_Coupling_Geometry:
            return "Coupling_Geometry";
        case GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Curve_On_Surface_Geometry:
            return "Quadrature_Point_Curve_On_Surface_Geometry";
        case GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Surface_In_Volume_Geometry:
            return "Quadrature_Point_Surface_In_Volume_Geometry";
        case GeometryData::KratosGeometryType::NumberOfGeometryTypes:
            return "NumberOfGeometryTypes";
        default:
            KRATOS_ERROR << "Geometry type not supported: " << static_cast<int>(TypeOfGeometry) << std::endl;
    };

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

double GeometryUtils::PointDistanceToLineSegment3D(
    const Point& rLinePoint1,
    const Point& rLinePoint2,
    const Point& rToPoint
    )
{
    const double epsilon = 1e-15; //1.0e-9;

    const array_1d<double,3> v1 = rLinePoint2 - rLinePoint1;
    const array_1d<double,3> v2 = rLinePoint1 - rToPoint;
    array_1d<double,3> v3;

    const double square_distance = inner_prod(v1,v1);

    if(square_distance < epsilon) // near zero length line
        return norm_2(v2); // we return the distance to the first point of line

    const double t = - inner_prod(v1,v2) / square_distance;

    if(t < 0.0) { // it is before point 1
        // We return the distance to point 1
        noalias(v3) = rLinePoint1 - rToPoint;

        return norm_2(v3);
    }

    if(t > 1.0) { // it is after point 2
        // We return the distance to point 2
        noalias(v3) = rLinePoint2 - rToPoint;

        return norm_2(v3);
    }

    // The projection point is between point 1 and 2 of the line segment
    noalias(v3) = rLinePoint1 * (1.0 - t) + rLinePoint2 * t;

    return norm_2(v3 - rToPoint);

}

/***********************************************************************************/
/***********************************************************************************/

double GeometryUtils::PointDistanceToTriangle3D(
    const Point& rTrianglePoint1,
    const Point& rTrianglePoint2,
    const Point& rTrianglePoint3,
    const Point& rPoint
    )
{
    const array_1d<double, 3> e0 = rTrianglePoint2 - rTrianglePoint1;
    const array_1d<double, 3> e1 = rTrianglePoint3 - rTrianglePoint1;
    const array_1d<double, 3> dd = rTrianglePoint1 - rPoint;

    const double a = inner_prod(e0, e0);
    const double b = inner_prod(e0, e1);
    const double c = inner_prod(e1, e1);
    const double d = inner_prod(e0, dd);
    const double e = inner_prod(e1, dd);
    const double f = inner_prod(dd, dd);

    const double det = a*c-b*b;
    double s = b*e-c*d;
    double t = b*d-a*e;

    double square_distance = 0.0;

    if ( s + t <= det ) {
        if ( s < 0.0 ) {
            if ( t < 0.0 ) { // region 4
                if (d < 0) {
                    t = 0;
                    if (-d >= a) {
                        s = 1;
                        square_distance = a + 2*d + f;
                    } else {
                        s = -d/a;
                        square_distance = d*s + f;
                    }
                } else {
                    s = 0;
                    if (e >= 0) {
                        t = 0;
                        square_distance = f;
                    } else {
                        if (-e >= c) {
                            t = 1;
                            square_distance = c + 2*e + f;
                        } else {
                            t = -e/c;
                            square_distance = e*t + f;
                        }
                    }
                }
            } else { // region 3
                s = 0.0;
                if(e >= 0.0) {
                    t = 0.0;
                    square_distance = f;
                } else {
                    if (-e >= c) {
                        t = 1.00;
                        square_distance = c + 2*e +f;
                    } else {
                        t = -e/c;
                        square_distance = e*t + f;
                    }
                }

            }
        } else if ( t < 0.00 ) { // region 5
            t = 0;
            if (d >= 0) {
                s = 0;
                square_distance = f;
            } else {
                if (-d >= a) {
                    s = 1;
                    square_distance = a + 2.0 * d + f;
                } else {
                    s = -d / a;
                    square_distance = d * s + f;
                }
            }
        } else { // region 0
            const double inv_det = 1.0 / det;
            s *= inv_det;
            t *= inv_det;
            square_distance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
        }
    } else {
        if ( s < 0.00 ) {
            // Region 2
            const double temp0 = b + d;
            const double temp1 = c + e;
            if (temp1 > temp0)  { // Minimum on edge s+t=1
                const double numer = temp1 - temp0;
                const double denom = a - 2*b + c;
                if(numer >= denom) {
                    s = 1.0;
                    t = 0.0;
                    square_distance = a + 2*d + f;
                } else {
                    s = numer/denom;
                    t = 1.0 - s;
                    square_distance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
                }
            } else { // Minimum on edge s=0
                s = 0.0;
                if(temp1 <= 0.0) {
                    t = 1;
                    square_distance = c + 2*e + f;
                } else {
                    if(e >= 0.0) {
                        t = 0.0;
                        square_distance = f;
                    } else {
                        t = -e/c;
                        square_distance = e*t + f;
                    }
                }
            }
        } else if ( t < 0.0 ) {
            // Region 6
            const double temp0 = b + e;
            const double temp1 = a + d;
            if (temp1 > temp0) {
                const double numer = temp1 - temp0;
                const double denom = a - 2*b + c;
                if (numer >= denom) {
                    s = 0.0;
                    t = 1.0;
                    square_distance = c + 2*e + f;
                } else {
                    t = numer/denom;
                    s = 1.0 - t;
                    square_distance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
                }
            } else {
                t = 0.0;
                if (temp1 <= 0.0) {
                    s = 1;
                    square_distance = a + 2*d + f;
                } else {
                    if(d >= 0.0) {
                        s = 0.0;
                        square_distance = f;
                    } else {
                        s = -d/a;
                        square_distance = d*s + f;
                    }
                }
            }
        } else {
            // Region 1
            const double numer = c + e - b - d;

            if (numer <= 0.0) {
                s = 0.0;
                t = 1.0;
                square_distance = c + 2.0 * e + f;
            } else {
                const double denom = a - 2.0 * b + c;
                if (numer >= denom) {
                    s = 1.0;
                    t = 0.0;
                    square_distance = a + 2.0 * d + f;
                } else {
                    s = numer / denom;
                    t = 1.0 - s;
                    square_distance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
                }
            }
        }
    }

    if(square_distance < 0.0)
        return 0.0; // avoiding -0 case!!

    return std::sqrt(square_distance);
}

/***********************************************************************************/
/***********************************************************************************/

double GeometryUtils::PointDistanceToTriangle3D(
    const Point& rTrianglePoint1,
    const Point& rTrianglePoint2,
    const Point& rTrianglePoint3,
    const Point& rTrianglePoint4,
    const Point& rTrianglePoint5,
    const Point& rTrianglePoint6,
    const Point& rPoint
    )
{
    std::array<double, 4> distances;
    distances[0] = GeometryUtils::PointDistanceToTriangle3D(rTrianglePoint1, rTrianglePoint4, rTrianglePoint6, rPoint);
    distances[1] = GeometryUtils::PointDistanceToTriangle3D(rTrianglePoint4, rTrianglePoint2, rTrianglePoint5, rPoint);
    distances[2] = GeometryUtils::PointDistanceToTriangle3D(rTrianglePoint6, rTrianglePoint5, rTrianglePoint3, rPoint);
    distances[3] = GeometryUtils::PointDistanceToTriangle3D(rTrianglePoint4, rTrianglePoint5, rTrianglePoint6, rPoint);
    const auto min = std::min_element(distances.begin(), distances.end());
    return *min;
}

/***********************************************************************************/
/***********************************************************************************/

double GeometryUtils::PointDistanceToQuadrilateral3D(
    const Point& rQuadrilateralPoint1,
    const Point& rQuadrilateralPoint2,
    const Point& rQuadrilateralPoint3,
    const Point& rQuadrilateralPoint4,
    const Point& rPoint
    )
{
    const double distance_1 = GeometryUtils::PointDistanceToTriangle3D(rQuadrilateralPoint1, rQuadrilateralPoint2, rQuadrilateralPoint3, rPoint);
    const double distance_2 = GeometryUtils::PointDistanceToTriangle3D(rQuadrilateralPoint3, rQuadrilateralPoint4, rQuadrilateralPoint1, rPoint);
    return std::min(distance_1, distance_2);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TDataType>
void GeometryUtils::EvaluateHistoricalVariableValueAtGaussPoint(
    TDataType& rOutput,
    const GeometryType& rGeometry,
    const Variable<TDataType>& rVariable,
    const Vector& rGaussPointShapeFunctionValues,
    const int Step)
{
    KRATOS_TRY

    const SizeType number_of_nodes = rGeometry.PointsNumber();

    noalias(rOutput) = rGeometry[0].FastGetSolutionStepValue(rVariable, Step) *
                       rGaussPointShapeFunctionValues[0];

    for (SizeType i_node = 1; i_node < number_of_nodes; ++i_node)
    {
        noalias(rOutput) += rGeometry[i_node].FastGetSolutionStepValue(rVariable, Step) *
                            rGaussPointShapeFunctionValues[i_node];
    }

    KRATOS_CATCH("");
}

template <>
void KRATOS_API(KRATOS_CORE) GeometryUtils::EvaluateHistoricalVariableValueAtGaussPoint<double>(
    double& rOutput,
    const GeometryType& rGeometry,
    const Variable<double>& rVariable,
    const Vector& rGaussPointShapeFunctionValues,
    const int Step)
{
    KRATOS_TRY

    const SizeType number_of_nodes = rGeometry.PointsNumber();

    rOutput = rGeometry[0].FastGetSolutionStepValue(rVariable, Step) *
              rGaussPointShapeFunctionValues[0];
    for (SizeType i_node = 1; i_node < number_of_nodes; ++i_node)
    {
        rOutput += rGeometry[i_node].FastGetSolutionStepValue(rVariable, Step) *
                   rGaussPointShapeFunctionValues[i_node];
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void GeometryUtils::EvaluateHistoricalVariableGradientAtGaussPoint(
    array_1d<double, 3>& rOutput,
    const GeometryType& rGeometry,
    const Variable<double>& rVariable,
    const Matrix& rGaussPointShapeFunctionDerivativeValues,
    const int Step)
{
    noalias(rOutput) = ZeroVector(3);
    const SizeType number_of_nodes = rGeometry.PointsNumber();
    const SizeType dimension = rGaussPointShapeFunctionDerivativeValues.size2();

    for (SizeType a = 0; a < number_of_nodes; ++a)
    {
        const double value = rGeometry[a].FastGetSolutionStepValue(rVariable, Step);
        for (SizeType i = 0; i < dimension; ++i)
            rOutput[i] += rGaussPointShapeFunctionDerivativeValues(a, i) * value;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void GeometryUtils::EvaluateHistoricalVariableGradientAtGaussPoint(
    BoundedMatrix<double, 3, 3>& rOutput,
    const GeometryType& rGeometry,
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rGaussPointShapeFunctionDerivativeValues,
    const int Step)
{
    noalias(rOutput) = ZeroMatrix(3, 3);
    const SizeType number_of_nodes = rGeometry.PointsNumber();
    const SizeType dimension = rGaussPointShapeFunctionDerivativeValues.size2();

    for (SizeType a = 0; a < number_of_nodes; ++a)
    {
        const array_1d<double, 3>& r_value =
            rGeometry[a].FastGetSolutionStepValue(rVariable, Step);
        for (SizeType i = 0; i < dimension; ++i)
        {
            for (SizeType j = 0; j < dimension; ++j)
            {
                rOutput(i, j) +=
                    rGaussPointShapeFunctionDerivativeValues(a, j) * r_value[i];
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryUtils::ProjectedIsInside(
    const GeometryType& rGeometry,
    const GeometryType::CoordinatesArrayType& rPointGlobalCoordinates,
    GeometryType::CoordinatesArrayType& rResult,
    const double Tolerance
    )
{
    // We compute the distance, if it is not in the pane we
    const Point point_to_project(rPointGlobalCoordinates);
    Point point_projected;
    double distance = 0.0;
    if (rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Line2D2) {
        distance = GeometricalProjectionUtilities::FastProjectOnLine2D(rGeometry, point_to_project, point_projected);
    } else if (rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3) {
        // We compute the normal to check the normal distances between the point and the triangles, so we can discard that is on the triangle
        const auto center = rGeometry.Center();
        const array_1d<double, 3> normal = rGeometry.UnitNormal(center);

        point_projected = GeometricalProjectionUtilities::FastProject(center, point_to_project, normal, distance);
    }

    // We check if we are on the plane
    if (std::abs(distance) > std::numeric_limits<double>::epsilon()) {
        if (std::abs(distance) > 1.0e-6 * rGeometry.Length()) {
            return false;
        }
    }

    return rGeometry.IsInside(rPointGlobalCoordinates, rResult, Tolerance);
}

// Template instantiations
template void KRATOS_API(KRATOS_CORE) GeometryUtils::EvaluateHistoricalVariableValueAtGaussPoint<array_1d<double, 3>>(
    array_1d<double, 3>& rOutput,
    const GeometryType&,
    const Variable<array_1d<double, 3>>&,
    const Vector&,
    const int);

} // namespace Kratos.
