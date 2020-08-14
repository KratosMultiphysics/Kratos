//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_RANS_APPLICATION_CALCULATION_UTILITIES_H_INCLUDED)
#define KRATOS_RANS_APPLICATION_CALCULATION_UTILITIES_H_INCLUDED

// System includes
#include <cmath>
#include <tuple>

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
using ElementType = ModelPart::ElementType;
using ConditionType = ModelPart::ConditionType;
/// Geometry type (using with given NodeType)
using GeometryType = Geometry<NodeType>;

template<class TDataType>
using RefVariablePair = std::tuple<TDataType&, const Variable<TDataType>&>;

inline long double SoftMax(
    const long double value_1,
    const long double value_2)
{
    return std::max(value_1, value_2);
}

inline long double SoftPositive(
    const long double value)
{
    return SoftMax(value, 0.0);
}

void CalculateGeometryData(
    const GeometryType& rGeometry,
    const GeometryData::IntegrationMethod& rIntegrationMethod,
    Vector& rGaussWeights,
    Matrix& rNContainer,
    GeometryType::ShapeFunctionsGradientsType& rDN_DX);

void CalculateConditionGeometryData(
    const GeometryType& rGeometry,
    const GeometryData::IntegrationMethod& rIntegrationMethod,
    Vector& rGaussWeights,
    Matrix& rNContainer);

GeometryType::ShapeFunctionsGradientsType CalculateGeometryParameterDerivatives(
    const GeometryType& rGeometry,
    const GeometryData::IntegrationMethod& rIntegrationMethod);

template <std::size_t TDim>
void CalculateGeometryParameterDerivativesShapeSensitivity(
    BoundedMatrix<double, TDim, TDim>& rOutput,
    const ShapeParameter& rShapeDerivative,
    const Matrix& rDnDe,
    const Matrix& rDeDx);

template<class TDataType>
void UpdateValue(TDataType& rOutput, const TDataType& rInput);

template<class TDataType>
constexpr RefVariablePair<TDataType> inline VariableValuePairTie(
    TDataType& rValue,
    const Variable<TDataType>& rVariable)
{
    // need this to make sure rVariables are ties with const references since
    // std::tie only support reference variable tieing.
    return std::tie(rValue, rVariable);
}

template<class... TDataTypeArgs>
void inline InitializeVariablePair(
    RefVariablePair<TDataTypeArgs>&&... rValueVariablePairs,
    const double ShapeFunctionValue,
    const NodeType& rNode,
    const int Step
)
{
    int dummy[sizeof...(TDataTypeArgs)] = {(
        std::get<0>(rValueVariablePairs) =
            rNode.FastGetSolutionStepValue(std::get<1>(rValueVariablePairs), Step) * ShapeFunctionValue,
        0)...};
    // following line is used to ignore warning of unused_variable can be removed by using fold expressions in c++17
    *dummy = 0;
}

template<class... TDataTypeArgs>
void inline UpdateVariablePair(
    RefVariablePair<TDataTypeArgs>&&... rValueVariablePairs,
    const double ShapeFunctionValue,
    const NodeType& rNode,
    const int Step
)
{
    int dummy[sizeof...(TDataTypeArgs)] = {(
        UpdateValue<TDataTypeArgs>(
            std::get<0>(rValueVariablePairs),
            rNode.FastGetSolutionStepValue(std::get<1>(rValueVariablePairs), Step) * ShapeFunctionValue),
        0)...};
    // following line is used to ignore warning of unused_variable can be removed by using fold expressions in c++17
    *dummy = 0;
}

template <class... TDataTypeArgs>
void EvaluateInPoint(
    const GeometryType& rGeometry,
    const Vector& rShapeFunction,
    const int Step,
    RefVariablePair<TDataTypeArgs>&&... rValueVariablePairs)
{
    const int number_of_nodes = rGeometry.PointsNumber();

    InitializeVariablePair<TDataTypeArgs...>(
        std::forward<RefVariablePair<TDataTypeArgs>>(rValueVariablePairs)...,
        rShapeFunction[0], rGeometry[0], Step);
    for (int c = 1; c < number_of_nodes; ++c) {
        UpdateVariablePair<TDataTypeArgs...>(
            std::forward<RefVariablePair<TDataTypeArgs>>(rValueVariablePairs)...,
            rShapeFunction[c], rGeometry[c], Step);
    }
}

template <class... TDataTypeArgs>
void EvaluateInPoint(
    const GeometryType& rGeometry,
    const Vector& rShapeFunction,
    RefVariablePair<TDataTypeArgs>&&... rValueVariablePairs)
{
    EvaluateInPoint<TDataTypeArgs...>(
        rGeometry, rShapeFunction, 0,
        std::forward<RefVariablePair<TDataTypeArgs>>(rValueVariablePairs)...);
}

template <unsigned int TDim>
double CalculateMatrixTrace(
    const BoundedMatrix<double, TDim, TDim>& rMatrix);

template <unsigned int TDim>
void CalculateGradient(
    BoundedMatrix<double, TDim, TDim>& rOutput,
    const Geometry<ModelPart::NodeType>& rGeometry,
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rShapeDerivatives,
    const int Step = 0);

void CalculateGradient(
    array_1d<double, 3>& rOutput,
    const Geometry<ModelPart::NodeType>& rGeometry,
    const Variable<double>& rVariable,
    const Matrix& rShapeDerivatives,
    const int Step = 0);

double GetDivergence(
    const Geometry<ModelPart::NodeType>& rGeometry,
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rShapeDerivatives,
    const int Step = 0);

template <unsigned int TNumNodes>
void CalculateGaussSensitivities(
    BoundedVector<double, TNumNodes>& rGaussSensitivities,
    const BoundedVector<double, TNumNodes>& rNodalSensitivities,
    const Vector& rGaussShapeFunctions);

template <unsigned int TDim>
Vector GetVector(
    const array_1d<double, 3>& rVector);

Vector GetVector(
    const array_1d<double, 3>& rVector,
    const unsigned int Dim);

double KRATOS_API(RANS_APPLICATION) CalculateLogarithmicYPlusLimit(
    const double Kappa,
    const double Beta,
    const int MaxIterations = 20,
    const double Tolerance = 1e-6);

void CalculateYPlusAndUtau(
    double& rYPlus,
    double& rUTau,
    const double WallVelocity,
    const double WallHeight,
    const double KinematicViscosity,
    const double Kappa,
    const double Beta,
    const int MaxIterations = 20,
    const double Tolerance = 1e-6);

double CalculateWallHeight(
    const ConditionType& rCondition,
    const array_1d<double, 3>& rNormal);

array_1d<double, 3> CalculateWallVelocity(
    const ConditionType& rCondition);

bool IsWallFunctionActive(
    const ConditionType& rCondition);

bool IsInlet(
    const ConditionType& rCondition);

} // namespace RansCalculationUtilities

///@}

} // namespace Kratos

#endif // KRATOS_RANS_APPLICATION_CALCULATION_UTILITIES_H_INCLUDED defined