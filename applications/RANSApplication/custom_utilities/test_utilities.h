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
typedef ModelPart::ConditionType ConditionType;
typedef Geometry<NodeType> GeometryType;
typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;
typedef std::size_t IndexType;

bool IsNear(double ValueA, double ValueB, double RelTol = 1e-09, double AbsTol = 1e-12);

void CheckNear(double ValueA, double ValueB, double RelTol = 1e-09, double AbsTol = 1e-12);

void CheckNear(const Matrix& rA, const Matrix& rB, double RelTol = 1e-09, double AbsTol = 1e-12);

void CheckNear(const Vector& rA, const Vector& rB, double RelTol = 1e-09, double AbsTol = 1e-12);

template <class TClassType>
void CalculateResidual(Vector& residual, TClassType& rClassTypeObject, ProcessInfo& rProcessInfo);

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

template <typename TContainer>
void RunResidualSensitivityTest(
    ModelPart& rPrimalModelPart,
    ModelPart& rAdjointModelPart,
    std::vector<Process*>& rPrimalProcesses,
    std::vector<Process*>& rAdjointProcesses,
    std::function<void(ModelPart&)> UpdateVariablesInModelPart,
    std::function<void(Matrix&, typename TContainer::data_type&, ProcessInfo&)> CalculateElementResidualScalarSensitivity,
    std::function<double&(NodeType&)> PerturbVariable,
    const double Delta,
    const double Tolerance,
    const int DerivativesOffset = 0,
    const int EquationOffset = 0);

template <typename TContainer>
void RunResidualSensitivityTest(
    ModelPart& rPrimalModelPart,
    ModelPart& rAdjointModelPart,
    std::vector<Process*>& rPrimalProcesses,
    std::vector<Process*>& rAdjointProcesses,
    std::function<void(ModelPart&)> UpdateVariablesInModelPart,
    std::function<void(Matrix&, typename TContainer::data_type&, ProcessInfo&)> CalculateElementResidualVectorSensitivity,
    std::function<double&(NodeType&, const int)> PerturbVariable,
    const double Delta,
    const double Tolerance,
    const int DerivativesOffset = 0,
    const int EquationOffset = 0);

template <typename TContainer>
TContainer& GetContainerItems(ModelPart& rModelPart);

std::function<double&(NodeType&)> GetPerturbationMethod(const Variable<double>&);

std::function<double&(NodeType&, const int)> GetPerturbationMethod(
    const Variable<array_1d<double, 3>>&);

} // namespace RansModellingApplicationTestUtilities
} // namespace Kratos

#endif