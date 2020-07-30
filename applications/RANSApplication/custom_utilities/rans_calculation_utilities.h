//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_RANS_APPLICATION_CALCULATION_UTILITIES_H_INCLUDED)
#define KRATOS_RANS_APPLICATION_CALCULATION_UTILITIES_H_INCLUDED

// System includes
#include <cmath>

// Project includes
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/model_part.h"
#include "utilities/geometrical_sensitivity_utility.h"

namespace Kratos
{
///@name Kratos Globals
///@{

namespace RansCalculationUtilities
{
/// Node type
using NodeType = ModelPart::NodeType;

/// Geometry type (using with given NodeType)
using GeometryType = Geometry<NodeType>;

inline long double SoftMax(const long double value_1, const long double value_2)
{
    return std::max(value_1, value_2);
}

inline long double SoftPositive(const long double value)
{
    return SoftMax(value, 0.0);
}

void CalculateGeometryData(const GeometryType& rGeometry,
                           const GeometryData::IntegrationMethod& rIntegrationMethod,
                           Vector& rGaussWeights,
                           Matrix& rNContainer,
                           GeometryType::ShapeFunctionsGradientsType& rDN_DX);

GeometryType::ShapeFunctionsGradientsType CalculateGeometryParameterDerivatives(
    const GeometryType& rGeometry, const GeometryData::IntegrationMethod& rIntegrationMethod);

template <std::size_t TDim>
void CalculateGeometryParameterDerivativesShapeSensitivity(BoundedMatrix<double, TDim, TDim>& rOutput,
                                                           const ShapeParameter& rShapeDerivative,
                                                           const Matrix& rDnDe,
                                                           const Matrix& rDeDx);

double EvaluateInPoint(const GeometryType& rGeometry,
                       const Variable<double>& rVariable,
                       const Vector& rShapeFunction,
                       const int Step = 0);

array_1d<double, 3> EvaluateInPoint(const GeometryType& rGeometry,
                                    const Variable<array_1d<double, 3>>& rVariable,
                                    const Vector& rShapeFunction,
                                    const int Step = 0);

template <unsigned int TDim>
double CalculateMatrixTrace(const BoundedMatrix<double, TDim, TDim>& rMatrix);

template <unsigned int TDim>
void CalculateGradient(BoundedMatrix<double, TDim, TDim>& rOutput,
                       const Geometry<ModelPart::NodeType>& rGeometry,
                       const Variable<array_1d<double, 3>>& rVariable,
                       const Matrix& rShapeDerivatives,
                       const int Step = 0);

void CalculateGradient(array_1d<double, 3>& rOutput,
                       const Geometry<ModelPart::NodeType>& rGeometry,
                       const Variable<double>& rVariable,
                       const Matrix& rShapeDerivatives,
                       const int Step = 0);

template <unsigned int TDim>
Vector GetVector(const array_1d<double, 3>& rVector);

Vector GetVector(const array_1d<double, 3>& rVector, const unsigned int Dim);

double KRATOS_API(RANS_APPLICATION)
    CalculateLogarithmicYPlusLimit(const double Kappa,
                                   const double Beta,
                                   const int MaxIterations = 20,
                                   const double Tolerance = 1e-6);

} // namespace RansCalculationUtilities

///@}

} // namespace Kratos

#endif // KRATOS_RANS_APPLICATION_CALCULATION_UTILITIES_H_INCLUDED defined