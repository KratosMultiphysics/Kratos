//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//
//

#if !defined(KRATOS_RANS_TEST_UTILITIES_H_INCLUDED)
#define KRATOS_RANS_TEST_UTILITIES_H_INCLUDED

// System includes
#include <functional>

// External includes

// Project includes
#include "includes/model_part.h"
#include "processes/process.h"

namespace Kratos
{
namespace RansModellingApplicationTestUtilities
{
typedef ModelPart::NodeType NodeType;
typedef ModelPart::ElementType ElementType;
typedef Geometry<NodeType> GeometryType;
typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

void IsValuesRelativelyNear(const double ValueA, const double ValueB, const double Tolerance);

void IsMatricesSame(const Matrix& rA, const Matrix& rB, const double Tolerance);

void IsVectorsSame(const Vector& rA, const Vector& rB, const double Tolerance);

void CalculateResidual(Vector& residual, Element& rElement, ProcessInfo& rProcessInfo);

void GetElementData(Vector& rGaussWeights,
                    Matrix& rShapeFunctions,
                    ShapeFunctionDerivativesArrayType& rShapeFunctionDerivatives,
                    const ElementType& rElement);

void InitializeVariableWithValues(ModelPart& rModelPart,
                                  const Variable<double>& rVariable,
                                  const double Value,
                                  const std::size_t TimeStep = 0);

void InitializeVariableWithRandomValues(ModelPart& rModelPart,
                                        const Variable<double>& rVariable,
                                        const double MinValue,
                                        const double MaxValue,
                                        const std::size_t TimeSteps);

void InitializeVariableWithRandomValues(ModelPart& rModelPart,
                                        const Variable<array_1d<double, 3>>& rVariable,
                                        const double MinValue,
                                        const double MaxValue,
                                        const std::size_t TimeSteps);

void RunGaussPointScalarSensitivityTest(
    ModelPart& rModelPart,
    Process& rYPlusProcess,
    Process& rNutProcess,
    std::function<void(std::vector<double>&, const ElementType&, const Vector&, const Matrix&, const ProcessInfo&)> CalculatePrimalQuantities,
    std::function<void(std::vector<Vector>&, const ElementType&, const Vector&, const Matrix&, const ProcessInfo&)> CalculateSensitivities,
    std::function<void(ModelPart&)> UpdateVariablesInModelPart,
    std::function<double&(NodeType&)> PerturbVariable,
    const double Delta,
    const double Tolerance);

void RunGaussPointVectorSensitivityTest(
    ModelPart& rModelPart,
    Process& rYPlusProcess,
    Process& rNutProcess,
    std::function<void(std::vector<double>&, const ElementType&, const Vector&, const Matrix&, const ProcessInfo&)> CalculatePrimalQuantities,
    std::function<void(std::vector<Matrix>&, const ElementType&, const Vector&, const Matrix&, const ProcessInfo&)> CalculateSensitivities,
    std::function<void(ModelPart&)> UpdateVariablesInModelPart,
    std::function<double&(NodeType&, const int Dim)> PerturbVariable,
    const double Delta,
    const double Tolerance);

void RunElementResidualScalarSensitivityTest(
    ModelPart& rPrimalModelPart,
    ModelPart& rAdjointModelPart,
    Process& rPrimalYPlusProcess,
    Process& rPrimalNutProcess,
    Process& rAdjointYPlusProcess,
    Process& rAdjointNutProcess,
    Process& rYPlusSensitivitiesProcess,
    Process& rNutSensitivitiesProcess,
    std::function<void(ModelPart&)> UpdateVariablesInModelPart,
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> CalculateElementResidualScalarSensitivity,
    std::function<double&(NodeType&)> PerturbVariable,
    const double Delta,
    const double Tolerance,
    const int DerivativesOffset = 0,
    const int EquationOffset = 0);

void RunElementResidualVectorSensitivityTest(
    ModelPart& rPrimalModelPart,
    ModelPart& rAdjointModelPart,
    Process& rPrimalYPlusProcess,
    Process& rPrimalNutProcess,
    Process& rAdjointYPlusProcess,
    Process& rAdjointNutProcess,
    Process& rYPlusSensitivitiesProcess,
    Process& rNutSensitivitiesProcess,
    std::function<void(ModelPart&)> UpdateVariablesInModelPart,
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> CalculateElementResidualVectorSensitivity,
    std::function<double&(NodeType&, const int)> PerturbVariable,
    const double Delta,
    const double Tolerance,
    const int DerivativesOffset = 0,
    const int EquationOffset = 0);

void RunNodalScalarSensitivityTest(
    ModelPart& rModelPart,
    Process& rYPlusProcess,
    std::function<void(std::vector<double>&, const NodeType&, const ProcessInfo&)> CalculatePrimalQuantities,
    std::function<void(std::vector<Vector>&, const ElementType&, const ProcessInfo&)> CalculateSensitivities,
    std::function<void(ModelPart&)> UpdateVariablesInModelPart,
    std::function<double&(NodeType&)> PerturbVariable,
    const double Delta,
    const double Tolerance);

void RunNodalVectorSensitivityTest(
    ModelPart& rModelPart,
    Process& rYPlusProcess,
    std::function<void(std::vector<double>&, const NodeType&, const ProcessInfo&)> CalculatePrimalQuantities,
    std::function<void(std::vector<Matrix>&, const ElementType&, const ProcessInfo&)> CalculateSensitivities,
    std::function<void(ModelPart&)> UpdateVariablesInModelPart,
    std::function<double&(NodeType&, const int Dim)> PerturbVariable,
    const double Delta,
    const double Tolerance);
} // namespace RansModellingApplicationTestUtilities
} // namespace Kratos

#endif