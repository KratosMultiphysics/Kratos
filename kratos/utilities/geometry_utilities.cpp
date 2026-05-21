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
#include "utilities/math_utils.h"
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
    // Compute the differences component-wise for e0
    const double e0_0 = rTrianglePoint2[0] - rTrianglePoint1[0];
    const double e0_1 = rTrianglePoint2[1] - rTrianglePoint1[1];
    const double e0_2 = rTrianglePoint2[2] - rTrianglePoint1[2];

    // Compute the differences component-wise for e1
    const double e1_0 = rTrianglePoint3[0] - rTrianglePoint1[0];
    const double e1_1 = rTrianglePoint3[1] - rTrianglePoint1[1];
    const double e1_2 = rTrianglePoint3[2] - rTrianglePoint1[2];

    // Compute the differences component-wise for dd
    const double dd_0 = rTrianglePoint1[0] - rPoint[0];
    const double dd_1 = rTrianglePoint1[1] - rPoint[1];
    const double dd_2 = rTrianglePoint1[2] - rPoint[2];

    // Compute the inner products explicitly
    const double a = e0_0 * e0_0 + e0_1 * e0_1 + e0_2 * e0_2;
    const double b = e0_0 * e1_0 + e0_1 * e1_1 + e0_2 * e1_2;
    const double c = e1_0 * e1_0 + e1_1 * e1_1 + e1_2 * e1_2;
    const double d = e0_0 * dd_0 + e0_1 * dd_1 + e0_2 * dd_2;
    const double e = e1_0 * dd_0 + e1_1 * dd_1 + e1_2 * dd_2;
    const double f = dd_0 * dd_0 + dd_1 * dd_1 + dd_2 * dd_2;

    // Compute the determinant
    const double det = a*c-b*b;
    double s = b*e-c*d;
    double t = b*d-a*e;

    // Initialize the square distance
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

void GeometryUtils::ShapeFunctionsSecondDerivativesTransformOnAllIntegrationPoints(
        DenseVector<DenseVector<Matrix>>& rResult,
        const GeometryType& rGeometry,
        const GeometryType::IntegrationMethod& rIntegrationMethod )
{

    const unsigned int integration_points_number = rGeometry.IntegrationPointsNumber( rIntegrationMethod );

    if ( integration_points_number == 0 )
            KRATOS_ERROR << "This integration method is not supported" << std::endl;

    rResult.resize(integration_points_number, false);

    //calculating the local gradients
     DenseVector<Matrix> DN_DX;
    rGeometry.ShapeFunctionsIntegrationPointsGradients( DN_DX, rIntegrationMethod );
    for (IndexType i = 0; i < integration_points_number; i++) {
        rResult[i].resize(rGeometry.PointsNumber(), false);
        for (IndexType j = 0; j < rGeometry.PointsNumber(); j++)
            rResult[i][j].resize(rGeometry.WorkingSpaceDimension(), rGeometry.WorkingSpaceDimension(), false);
        GeometryUtils::ShapeFunctionsSecondDerivativesTransformOnIntegrationPoint(DN_DX[i],rGeometry,rGeometry.IntegrationPoints( rIntegrationMethod )[i], rResult[i]);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void GeometryUtils::ShapeFunctionsSecondDerivativesTransformOnIntegrationPoint(
      const Matrix& DN_DX,
      const GeometryType& rGeometry,
      const GeometryType::CoordinatesArrayType& rLocalIntegrationPointCoordinates,
      DenseVector<Matrix>& rResult)
{
    KRATOS_ERROR_IF_NOT(rGeometry.WorkingSpaceDimension() == rGeometry.LocalSpaceDimension())
        << "\'ShapeFunctionsIntegrationPointsSecondDerivatives\' is not defined for current geometry type as second derivatives are only defined in the local space." << std::endl;

    GeometryType::ShapeFunctionsSecondDerivativesType DDN_DDe;
    rGeometry.ShapeFunctionsSecondDerivatives(DDN_DDe, rLocalIntegrationPointCoordinates);

    Matrix A, Ainv;
    double DetA;
    Matrix J(rGeometry.WorkingSpaceDimension(),rGeometry.LocalSpaceDimension());
    rGeometry.Jacobian(J,rLocalIntegrationPointCoordinates);

    DenseVector<Matrix> aux(rGeometry.PointsNumber());
    Vector rhs, result;
    if (rGeometry.WorkingSpaceDimension() == 2){
        rhs.resize(3, false);
        result.resize(3, false);
    }
    else if (rGeometry.WorkingSpaceDimension() == 3){
        rhs.resize(6, false);
        result.resize(6, false);
    }

    for (IndexType i = 0; i < rGeometry.PointsNumber(); i++) {
        aux[i].resize(rGeometry.WorkingSpaceDimension(), rGeometry.WorkingSpaceDimension(), false);
    }

    if (rGeometry.WorkingSpaceDimension() == 2){
            A.resize(3,3,false);
            Ainv.resize(3,3,false);

            A(0,0) = J(0,0) * J(0,0);
            A(0,1) = J(1,0) * J(1,0);
            A(0,2) = 2.0 * J(0,0) * J(1,0);

            A(1,0) = J(0,1) * J(0,1);
            A(1,1) = J(1,1) * J(1,1);
            A(1,2) = 2.0 * J(0,1) * J(1,1);

            A(2,0) = J(0,0) * J(0,1);
            A(2,1) = J(1,0) * J(1,1);
            A(2,2) = J(0,0) * J(1,1) + J(0,1) * J(1,0);
        }
        else if (rGeometry.WorkingSpaceDimension() == 3){
            A.resize(6,6,false);
            Ainv.resize(6,6,false);

            A(0,0) = J(0,0) * J(0,0);
            A(0,1) = J(1,0) * J(1,0);
            A(0,2) = J(2,0) * J(2,0);
            A(0,3) = 2.0 * J(0,0) * J(1,0);
            A(0,4) = 2.0 * J(1,0) * J(2,0);
            A(0,5) = 2.0 * J(0,0) * J(2,0);

            A(1,0) = J(0,1) * J(0,1);
            A(1,1) = J(1,1) * J(1,1);
            A(1,2) = J(2,1) * J(2,1);
            A(1,3) = 2.0 * J(0,1) * J(1,1);
            A(1,4) = 2.0 * J(1,1) * J(2,1);
            A(1,5) = 2.0 * J(0,1) * J(2,1);

            A(2,0) = J(0,2) * J(0,2);
            A(2,1) = J(1,2) * J(1,2);
            A(2,2) = J(2,2) * J(2,2);
            A(2,3) = 2.0 * J(0,2) * J(1,2);
            A(2,4) = 2.0 * J(1,2) * J(2,2);
            A(2,5) = 2.0 * J(0,2) * J(2,2);

            A(3,0) = J(0,0) * J(0,1);
            A(3,1) = J(1,0) * J(1,1);
            A(3,2) = J(2,0) * J(2,1);
            A(3,3) = J(0,0) * J(1,1) + J(0,1) * J(1,0);
            A(3,4) = J(1,0) * J(2,1) + J(1,1) * J(2,0);
            A(3,5) = J(0,0) * J(2,1) + J(0,1) * J(2,0);

            A(4,0) = J(0,1) * J(0,2);
            A(4,1) = J(1,1) * J(1,2);
            A(4,2) = J(2,1) * J(2,2);
            A(4,3) = J(0,1) * J(1,2) + J(0,2) * J(1,1);
            A(4,4) = J(1,1) * J(2,2) + J(1,2) * J(2,1);
            A(4,5) = J(0,1) * J(2,2) + J(0,2) * J(2,1);

            A(5,0) = J(0,0) * J(0,2);
            A(5,1) = J(1,0) * J(1,2);
            A(5,2) = J(2,0) * J(2,2);
            A(5,3) = J(0,0) * J(1,2) + J(0,2) * J(1,0);
            A(5,4) = J(1,0) * J(2,2) + J(1,2) * J(2,0);
            A(5,5) = J(0,0) * J(2,2) + J(0,2) * J(2,0);

        }
        MathUtils<double>::InvertMatrix( A, Ainv, DetA );
        DenseVector<Matrix> H(rGeometry.WorkingSpaceDimension());
        for (unsigned int d = 0; d < rGeometry.WorkingSpaceDimension(); ++d)
            H[d] = ZeroMatrix(rGeometry.LocalSpaceDimension(),rGeometry.LocalSpaceDimension());

        for (IndexType p = 0; p < rGeometry.PointsNumber(); ++p) {
            const array_1d<double, 3>& r_coordinates = rGeometry[p].Coordinates();
            for (unsigned int d=0 ; d < rGeometry.WorkingSpaceDimension(); ++d){
                for (unsigned int e=0 ; e < rGeometry.WorkingSpaceDimension(); ++e){
                    for (unsigned int f=0 ; f < rGeometry.WorkingSpaceDimension(); ++f){
                        H[d](e,f) += r_coordinates[d] * DDN_DDe[p](e,f);
                    }
                }
            }
        }
        for (IndexType p = 0; p < rGeometry.PointsNumber(); ++p) {
            if (rGeometry.WorkingSpaceDimension() == 2){
                rhs[0] = DDN_DDe[p](0,0) - DN_DX(p,0) * H[0](0,0) - DN_DX(p,1) * H[1](0,0);
                rhs[1] = DDN_DDe[p](1,1) - DN_DX(p,0) * H[0](1,1) - DN_DX(p,1) * H[1](1,1);
                rhs[2] = DDN_DDe[p](0,1) - DN_DX(p,0) * H[0](0,1) - DN_DX(p,1) * H[1](0,1);
            }
            else if (rGeometry.WorkingSpaceDimension() == 3){
                rhs[0] = DDN_DDe[p](0,0) - DN_DX(p,0) * H[0](0,0) - DN_DX(p,1) * H[1](0,0) - DN_DX(p,2) * H[2](0,0);
                rhs[1] = DDN_DDe[p](1,1) - DN_DX(p,0) * H[0](1,1) - DN_DX(p,1) * H[1](1,1) - DN_DX(p,2) * H[2](0,0);
                rhs[2] = DDN_DDe[p](2,2) - DN_DX(p,0) * H[0](2,2) - DN_DX(p,1) * H[1](2,2) - DN_DX(p,2) * H[2](2,2);
                rhs[3] = DDN_DDe[p](0,1) - DN_DX(p,0) * H[0](0,1) - DN_DX(p,1) * H[1](0,1) - DN_DX(p,2) * H[2](0,1);
                rhs[4] = DDN_DDe[p](1,2) - DN_DX(p,0) * H[0](1,2) - DN_DX(p,1) * H[1](1,2) - DN_DX(p,2) * H[2](1,2);
                rhs[5] = DDN_DDe[p](0,2) - DN_DX(p,0) * H[0](0,2) - DN_DX(p,1) * H[1](0,2) - DN_DX(p,2) * H[2](0,2);
            }

            aux[p].resize(rGeometry.WorkingSpaceDimension(), rGeometry.WorkingSpaceDimension(), false );

            noalias(result) = prod(Ainv, rhs);

            aux[p](0,0) = result[0];
            aux[p](1,1) = result[1];
            if (rGeometry.WorkingSpaceDimension() == 2){
                aux[p](0,1) = result[2];
                aux[p](1,0) = result[2];
            }
            else if (rGeometry.WorkingSpaceDimension() == 3){
                aux[p](2,2) = result[2];
                aux[p](0,1) = result[3];
                aux[p](1,0) = result[3];
                aux[p](0,2) = result[5];
                aux[p](2,0) = result[5];
                aux[p](2,1) = result[4];
                aux[p](1,2) = result[4];
            }
        }
        noalias(rResult) = aux;
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

/***********************************************************************************/
/***********************************************************************************/

bool GeometryUtils::TriangleBoxOverlap(
    const Point& rBoxCenter,
    const Point& rBoxHalfSize,
    const Point& rVertex0,
    const Point& rVertex1,
    const Point& rVertex2
    )
{
    double abs_ex, abs_ey, abs_ez, distance;
    array_1d<double,3 > vert0, vert1, vert2;
    array_1d<double,3 > edge0, edge1, edge2, normal;
    std::pair<double, double> min_max;

    // move everything so that the boxcenter is in (0,0,0)
    noalias(vert0) = rVertex0 - rBoxCenter;
    noalias(vert1) = rVertex1 - rBoxCenter;
    noalias(vert2) = rVertex2 - rBoxCenter;

    // compute triangle edges
    noalias(edge0) = vert1 - vert0;
    noalias(edge1) = vert2 - vert1;
    noalias(edge2) = vert0 - vert2;

    // Bullet 3:
    // test the 9 tests first (this was faster)
    abs_ex = std::abs(edge0[0]);
    abs_ey = std::abs(edge0[1]);
    abs_ez = std::abs(edge0[2]);
    if (AxisTestX(edge0[1],edge0[2],abs_ey,abs_ez,vert0,vert2,rBoxHalfSize)) return false;
    if (AxisTestY(edge0[0],edge0[2],abs_ex,abs_ez,vert0,vert2,rBoxHalfSize)) return false;
    if (AxisTestZ(edge0[0],edge0[1],abs_ex,abs_ey,vert0,vert2,rBoxHalfSize)) return false;

    abs_ex = std::abs(edge1[0]);
    abs_ey = std::abs(edge1[1]);
    abs_ez = std::abs(edge1[2]);
    if (AxisTestX(edge1[1],edge1[2],abs_ey,abs_ez,vert1,vert0,rBoxHalfSize)) return false;
    if (AxisTestY(edge1[0],edge1[2],abs_ex,abs_ez,vert1,vert0,rBoxHalfSize)) return false;
    if (AxisTestZ(edge1[0],edge1[1],abs_ex,abs_ey,vert1,vert0,rBoxHalfSize)) return false;

    abs_ex = std::abs(edge2[0]);
    abs_ey = std::abs(edge2[1]);
    abs_ez = std::abs(edge2[2]);
    if (AxisTestX(edge2[1],edge2[2],abs_ey,abs_ez,vert2,vert1,rBoxHalfSize)) return false;
    if (AxisTestY(edge2[0],edge2[2],abs_ex,abs_ez,vert2,vert1,rBoxHalfSize)) return false;
    if (AxisTestZ(edge2[0],edge2[1],abs_ex,abs_ey,vert2,vert1,rBoxHalfSize)) return false;

    // Bullet 1:
    //  first test overlap in the {x,y,z}-directions
    //  find min, max of the triangle for each direction, and test for
    //  overlap in that direction -- this is equivalent to testing a minimal
    //  AABB around the triangle against the AABB

    // test in X-direction
    min_max = std::minmax({vert0[0], vert1[0], vert2[0]});
    if(min_max.first>rBoxHalfSize[0] || min_max.second<-rBoxHalfSize[0]) return false;

    // test in Y-direction
    min_max = std::minmax({vert0[1], vert1[1], vert2[1]});
    if(min_max.first>rBoxHalfSize[1] || min_max.second<-rBoxHalfSize[1]) return false;

    // test in Z-direction
    min_max = std::minmax({vert0[2], vert1[2], vert2[2]});
    if(min_max.first>rBoxHalfSize[2] || min_max.second<-rBoxHalfSize[2]) return false;

    // Bullet 2:
    //  test if the box intersects the plane of the triangle
    //  compute plane equation of triangle: normal*x+distance=0
    MathUtils<double>::CrossProduct(normal, edge0, edge1);
    distance = -inner_prod(normal, vert0);
    if(!PlaneBoxOverlap(normal, distance, rBoxHalfSize)) return false;

    return true;  // box and triangle overlaps
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryUtils::AxisTestX(
    const double EdgeY,
    const double EdgeZ,
    const double AbsEdgeY,
    const double AbsEdgeZ,
    const array_1d<double,3>& rVertA,
    const array_1d<double,3>& rVertC,
    const Point& rBoxHalfSize
    )
{
    double proj_a, proj_c, rad;
    proj_a = EdgeY*rVertA[2] - EdgeZ*rVertA[1];
    proj_c = EdgeY*rVertC[2] - EdgeZ*rVertC[1];
    std::pair<double, double> min_max = std::minmax(proj_a, proj_c);

    rad = AbsEdgeZ*rBoxHalfSize[1] + AbsEdgeY*rBoxHalfSize[2];

    if(min_max.first>rad || min_max.second<-rad) return true;
    else return false;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryUtils::AxisTestY(
    const double EdgeX,
    const double EdgeY,
    const double AbsEdgeX,
    const double AbsEdgeZ,
    const array_1d<double,3>& rVertA,
    const array_1d<double,3>& rVertC,
    const Point& rBoxHalfSize
    )
{
    double proj_a, proj_c, rad;
    proj_a = EdgeY*rVertA[0] - EdgeX*rVertA[2];
    proj_c = EdgeY*rVertC[0] - EdgeX*rVertC[2];
    std::pair<double, double> min_max = std::minmax(proj_a, proj_c);

    rad = AbsEdgeZ*rBoxHalfSize[0] + AbsEdgeX*rBoxHalfSize[2];

    if(min_max.first>rad || min_max.second<-rad) return true;
    else return false;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryUtils::AxisTestZ(
    const double EdgeX,
    const double EdgeY,
    const double AbsEdgeX,
    const double AbsEdgeY,
    const array_1d<double,3>& rVertA,
    const array_1d<double,3>& rVertC,
    const Point& rBoxHalfSize
    )
{
    double proj_a, proj_c, rad;
    proj_a = EdgeX*rVertA[1] - EdgeY*rVertA[0];
    proj_c = EdgeX*rVertC[1] - EdgeY*rVertC[0];
    std::pair<double, double> min_max = std::minmax(proj_a, proj_c);

    rad = AbsEdgeY*rBoxHalfSize[0] + AbsEdgeX*rBoxHalfSize[1];

    if(min_max.first>rad || min_max.second<-rad) return true;
    else return false;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryUtils::PlaneBoxOverlap(
    const array_1d<double,3>& rNormal,
    const double Distance,
    const array_1d<double,3>& rMaxBox
    )
{
    array_1d<double,3> vmin, vmax;
    for(int q = 0; q < 3; q++) {
        if(rNormal[q] > 0.00) {
            vmin[q] = -rMaxBox[q];
            vmax[q] =  rMaxBox[q];
        } else {
            vmin[q] =  rMaxBox[q];
            vmax[q] = -rMaxBox[q];
        }
    }
    if(inner_prod(rNormal, vmin) + Distance >  0.0) return false;
    if(inner_prod(rNormal, vmax) + Distance >= 0.0) return true;

    return false;
}

} // namespace Kratos.