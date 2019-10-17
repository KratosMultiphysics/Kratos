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

void CreateModelPartElements(ModelPart& rModelPart, std::string ElementName);

void CreateEquationIds(ModelPart& rModelPart);

void InitializeNodalVariables(ModelPart& rModelPart);

void GenerateRansEvmKEpsilonElementTestModelPart(ModelPart& rModelPart, std::string ElementName);

void GenerateRansEvmKEpsilonConditionTestModelPart(ModelPart& rModelPart,
                                                   std::string ConditionName);

void UpdateVariablesInModelPart(ModelPart& rModelPart);

void UpdateVariablesInModelPartLowRe(ModelPart& rModelPart);

void CalculatePrimalQuantities(std::vector<double>& rValues,
                               const ElementType& rElement,
                               const Vector& rGaussShapeFunctions,
                               const Matrix& rGaussShapeFunctionDerivatives,
                               const ProcessInfo& rCurrentProcessInfo);

void ReadNodalDataFromElement(Vector& rYPlus,
                              Vector& rTKE,
                              Vector& rEpsilon,
                              Vector& rNut,
                              Vector& rFmu,
                              const Element& rElement);

void InitializeYPlus(ModelPart& rModelPart);

} // namespace RansEvmKEpsilonModel
} // namespace Testing
} // namespace Kratos
#endif
