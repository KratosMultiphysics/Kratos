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
// TODO: Remove this after std14 upgrade
// this is not required once we upgrade to std14 or better :)
namespace std14
{
template <std::size_t...>
struct index_sequence {
};

template <std::size_t N, std::size_t... Next>
struct indexSequenceHelper : public indexSequenceHelper<N - 1U, N - 1U, Next...> {
};

template <std::size_t... Next>
struct indexSequenceHelper<0U, Next...> {
    using type = index_sequence<Next...>;
};

template <std::size_t N>
using make_index_sequence = typename indexSequenceHelper<N>::type;
} // namespace std14

/// Node type
using NodeType = ModelPart::NodeType;
using ElementType = ModelPart::ElementType;
using ConditionType = ModelPart::ConditionType;
/// Geometry type (using with given NodeType)
using GeometryType = Geometry<NodeType>;

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

double EvaluateInPoint(
    const GeometryType& rGeometry,
    const Variable<double>& rVariable,
    const Vector& rShapeFunction,
    const int Step = 0);

template<class TDataType>
void UpdateValue(TDataType& rOutput, const TDataType& rInput);

template<class... TDataTypeArgs, std::size_t... Is>
void inline InitializePartialGaussPointValues(
    std::tuple<TDataTypeArgs&...> rValues,
    const double ShapeFunctionValue,
    const NodeType& rNode,
    const int Step,
    const std14::index_sequence<Is...>&,
    const Variable<TDataTypeArgs>&... rVariables)
{
    int dummy[sizeof...(TDataTypeArgs)] = {(
        std::get<Is>(rValues) = rNode.FastGetSolutionStepValue(rVariables, Step) * ShapeFunctionValue,
        0)...};
    // following line is used to ignore warning of unused_variable
    *dummy = 0;
}

template<class... TDataTypeArgs, std::size_t... Is>
void inline UpdatePartialGaussPointValues(
    std::tuple<TDataTypeArgs&...> rValues,
    const double ShapeFunctionValue,
    const NodeType& rNode,
    const int Step,
    const std14::index_sequence<Is...>&,
    const Variable<TDataTypeArgs>&... rVariables)
{
    int dummy[sizeof...(TDataTypeArgs)] = {(
        UpdateValue<TDataTypeArgs>(std::get<Is>(rValues),
                                   rNode.FastGetSolutionStepValue(rVariables, Step) * ShapeFunctionValue),
        0)...};
    // following line is used to ignore warning of unused_variable
    *dummy = 0;
}

template <class... TDataTypeArgs>
void EvaluateInPoint(
    std::tuple<TDataTypeArgs&...> rValues,
    const GeometryType& rGeometry,
    const Vector& rShapeFunction,
    const int Step,
    const Variable<TDataTypeArgs>&... rVariables)
{
    const int number_of_nodes = rGeometry.PointsNumber();
    const auto indexed_sequence = std14::make_index_sequence<sizeof...(TDataTypeArgs)>();

    InitializePartialGaussPointValues<TDataTypeArgs...>(
        rValues, rShapeFunction[0], rGeometry[0], Step, indexed_sequence, rVariables...);
    for (int c = 1; c < number_of_nodes; ++c) {
        UpdatePartialGaussPointValues<TDataTypeArgs...>(
            rValues, rShapeFunction[c], rGeometry[c], Step, indexed_sequence, rVariables...);
    }
}

template <class... TDataTypeArgs>
void EvaluateInPoint(
    std::tuple<TDataTypeArgs&...> rValues,
    const GeometryType& rGeometry,
    const Vector& rShapeFunction,
    const Variable<TDataTypeArgs>&... rVariables)
{
    EvaluateInPoint<TDataTypeArgs...>(rValues, rGeometry, rShapeFunction, 0, rVariables...);
}

array_1d<double, 3> EvaluateInPoint(
    const GeometryType& rGeometry,
    const Variable<array_1d<double, 3>>& rVariable,
    const Vector& rShapeFunction,
    const int Step = 0);

template <typename TDataType>
TDataType EvaluateInParentCenter(
    const Variable<TDataType>& rVariable,
    const ConditionType& rCondition,
    const int Step = 0);

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