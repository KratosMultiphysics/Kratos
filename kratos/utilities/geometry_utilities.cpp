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
            return "Triangle3D3";
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
            KRATOS_ERROR << "Geometry type not supported" << std::endl;
            return "NumberOfGeometryTypes";
        default:
            KRATOS_ERROR << "Geometry type not supported: " << static_cast<int>(TypeOfGeometry) << std::endl;
    };

    KRATOS_CATCH("");
}

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

// template instantiations

template void KRATOS_API(KRATOS_CORE) GeometryUtils::EvaluateHistoricalVariableValueAtGaussPoint<array_1d<double, 3>>(
    array_1d<double, 3>& rOutput,
    const GeometryType&,
    const Variable<array_1d<double, 3>>&,
    const Vector&,
    const int);

} // namespace Kratos.
