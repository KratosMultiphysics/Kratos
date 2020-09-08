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

// TODO: Move this to core GeometryUtils
void CalculateGeometryData(
    const GeometryType& rGeometry,
    const GeometryData::IntegrationMethod& rIntegrationMethod,
    Vector& rGaussWeights,
    Matrix& rNContainer,
    GeometryType::ShapeFunctionsGradientsType& rDN_DX);

// TODO: Move this to core GeometryUtils
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

template <class... TRefVariableValuePairArgs>
void EvaluateInPoint(
    const GeometryType& rGeometry,
    const Vector& rShapeFunction,
    const int Step,
    const TRefVariableValuePairArgs&... rValueVariablePairs)
{
    KRATOS_TRY

    const int number_of_nodes = rGeometry.PointsNumber();

    const auto& r_node = rGeometry[0];
    const double shape_function_value = rShapeFunction[0];

    int dummy[sizeof...(TRefVariableValuePairArgs)] = {(
        std::get<0>(rValueVariablePairs) =
            r_node.FastGetSolutionStepValue(std::get<1>(rValueVariablePairs), Step) * shape_function_value,
        0)...};

    // this can be removed with fold expressions in c++17
    *dummy = 0;

    for (int c = 1; c < number_of_nodes; ++c) {
        const auto& r_node = rGeometry[c];
        const double shape_function_value = rShapeFunction[c];

        int dummy[sizeof...(TRefVariableValuePairArgs)] = {(
            UpdateValue<typename std::remove_reference<typename std::tuple_element<0, TRefVariableValuePairArgs>::type>::type>(
                std::get<0>(rValueVariablePairs),
                r_node.FastGetSolutionStepValue(std::get<1>(rValueVariablePairs), Step) * shape_function_value),
            0)...};

        // this can be removed with fold expressions in c++17
        *dummy = 0;
    }

    KRATOS_CATCH("");
}

template <class... TRefVariableValuePairArgs>
void inline EvaluateInPoint(
    const GeometryType& rGeometry,
    const Vector& rShapeFunction,
    const TRefVariableValuePairArgs&... rValueVariablePairs)
{
    EvaluateInPoint<TRefVariableValuePairArgs...>(rGeometry, rShapeFunction, 0, rValueVariablePairs...);
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

template <class TContainerType>
TContainerType& GetContainer(ModelPart& rModelPart);

/**
 * @brief Calculates number of neighbours
 *
 * Calculates number of neighbours for a given entities in given model part.
 * The number of neighbours will be stored in rOutputVariable
 *
 * @tparam TContainerType       Entities type (ConditionsContainerType, ElementsContainerType)
 * @param rModelPart            Model part to look for neighbours
 * @param rOutputVariable       Variable to store number of neighbour entities.
 */
template<class TContainerType>
void CalculateNumberOfNeighbourEntities(
    ModelPart& rModelPart,
    const Variable<double>& rOutputVariable);

} // namespace RansCalculationUtilities

///@}

} // namespace Kratos

#endif // KRATOS_RANS_APPLICATION_CALCULATION_UTILITIES_H_INCLUDED defined