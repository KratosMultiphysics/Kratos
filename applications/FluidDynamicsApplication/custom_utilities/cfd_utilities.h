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

#if !defined(KRATOS_CFD_UTILITIES_H_INCLUDED)
#define KRATOS_CFD_UTILITIES_H_INCLUDED

// System includes
#include <iostream>
#include <string>
#include <tuple>

// External includes

// Project includes
#include "containers/array_1d.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{
namespace CFDUtilities
{
using NodeType = ModelPart::NodeType;
using ElementType = ModelPart::ElementType;
using ConditionType = ModelPart::ConditionType;
using GeometryType = Geometry<NodeType>;

void CalculateConditionGeometryData(const GeometryType& rGeometry,
                                    const GeometryData::IntegrationMethod& rIntegrationMethod,
                                    Vector& rGaussWeights,
                                    Matrix& rNContainer);

template <unsigned int TDim>
void CalculateConditionNormal(array_1d<double, 3>& rNormal, const ConditionType& rCondition);

double CalculateConditionWallHeight(const ConditionType& rCondition,
                                    const array_1d<double, 3>& rNormal);

template <typename TDataType>
TDataType EvaluateInPoint(const GeometryType& rGeometry,
                          const Variable<TDataType>& rVariable,
                          const Vector& rShapeFunction,
                          const int Step = 0)
{
    const unsigned int number_of_nodes = rGeometry.PointsNumber();
    TDataType value =
        rGeometry[0].FastGetSolutionStepValue(rVariable, Step) * rShapeFunction[0];
    for (unsigned int c = 1; c < number_of_nodes; ++c)
    {
        value += rGeometry[c].FastGetSolutionStepValue(rVariable, Step) *
                 rShapeFunction[c];
    }

    return value;
}

void KRATOS_API(FLUID_DYNAMICS_APPLICATION)
    CalculateNumberOfNeighbourConditions(ModelPart& rModelPart);

double KRATOS_API(FLUID_DYNAMICS_APPLICATION)
    CalculateLinearLogarithmicWallFunctionBasedYPlusLimit(const double VonKarman = 0.41,
                                                          const double WallSmoothness = 5.2,
                                                          const int MaxIterations = 20,
                                                          const double Tolerance = 1e-6);

double KRATOS_API(FLUID_DYNAMICS_APPLICATION)
    CalculateLinearLogarithmicWallFunctionBasedYPlusAndUtau(
        array_1d<double, 3>& rFrictionVelocity,
        const array_1d<double, 3>& rWallVelocity,
        const array_1d<double, 3>& rNormal,
        const double KinematicViscosity,
        const double WallHeight,
        const double VonKarman = 0.41,
        const double WallSmoothness = 5.2,
        const int MaxIterations = 20,
        const double Tolerance = 1e-6);

double KRATOS_API(FLUID_DYNAMICS_APPLICATION)
    CalculateReactionBasedYPlusUTau(array_1d<double, 3>& rFrictionVelocity,
                                    const array_1d<double, 3>& rReaction,
                                    const array_1d<double, 3>& rNormal,
                                    const double Density,
                                    const double KinematicViscosity,
                                    const double WallHeight);

void CalculateYPlusAndUTauForConditions(
    ModelPart& rModelPart,
    const Variable<double>& rKinematicViscosityVariable,
    const std::function<double(
        array_1d<double, 3>&, const GeometryType&, const array_1d<double, 3>&, const Vector&, const double, const double, const double)>&
        rYPlusAndUTauCalculationMethod);

void KRATOS_API(FLUID_DYNAMICS_APPLICATION) CalculateYPlusAndUTauForConditionsBasedOnReaction(
    ModelPart& rModelPart,
    const Variable<double>& rKinematicViscosityVariable,
    const Variable<array_1d<double, 3>>& rReactionVariable);

void KRATOS_API(FLUID_DYNAMICS_APPLICATION)
    CalculateYPlusAndUTauForConditionsBasedOnLinearLogarithmicWallFunction(
        ModelPart& rModelPart,
        const Variable<double>& rKinematicViscosityVariable,
        const double VonKarman = 0.41,
        const double WallSmoothness = 5.2,
        const int MaxIterations = 20,
        const double Tolerance = 1e-6);

template <typename TDataType>
void DistributeConditionVariableToNodes(ModelPart& rModelPart,
                                        const Variable<TDataType>& rVariable);

} // namespace CFDUtilities

} // namespace Kratos.

#endif // KRATOS_CFD_UTILITIES_H_INCLUDED  defined
