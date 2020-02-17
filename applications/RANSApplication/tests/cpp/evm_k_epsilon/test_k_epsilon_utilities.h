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

#if !defined(KRATOS_RANS_TEST_K_EPSILON_UTILITIES_H_INCLUDED)
#define KRATOS_RANS_TEST_K_EPSILON_UTILITIES_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/process_info.h"

// Application includes

namespace Kratos
{
namespace Testing
{
typedef ModelPart::NodeType NodeType;

typedef ModelPart::ElementType ElementType;

typedef ModelPart::ConditionType ConditionType;

typedef Geometry<NodeType> GeometryType;

typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

/**
 * Auxiliar function to generate a triangular element to be tested.
 */
namespace RansEvmKEpsilonModel
{
void AddVariablesToModelPart(ModelPart& rModelPart);

void InitializeProcessInfo(ModelPart& rModelPart);

void CreateModelPartNodes(ModelPart& rModelPart);

void CreateModelPartElements(ModelPart& rModelPart, const std::string& ElementName);

void CreateEquationIds(ModelPart& rModelPart);

void InitializeNodalVariables(ModelPart& rModelPart);

template <typename TContainerType>
void GenerateRansEvmKEpsilonTestModelPart(ModelPart& rModelPart,
                                          const std::string& TContainerDataTypeName);

void UpdateVariablesInModelPart(ModelPart& rModelPart);

void InitializeYPlus(ModelPart& rModelPart);

template <typename TDataType, typename TContainer>
void RunRansEvmKEpsilonTest(const std::string& PrimalName,
                            const std::string& AdjointName,
                            const Variable<TDataType>& rPerturbationVariable,
                            std::function<void(Matrix&, typename TContainer::data_type&, ProcessInfo&)> CalculateElementResidualScalarSensitivity,
                            const double Delta,
                            const double Tolerance,
                            const int DerivativesOffset = 0,
                            const int EquationOffset = 0);

} // namespace RansEvmKEpsilonModel
} // namespace Testing
} // namespace Kratos
#endif
