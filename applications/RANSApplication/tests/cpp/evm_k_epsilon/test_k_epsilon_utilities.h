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
#include "custom_processes/auxiliary_processes/rans_nut_k_epsilon_high_re_update_process.h"
#include "custom_utilities/test_utilities.h"

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

template <unsigned int TDim, unsigned int TNumNodes, typename TDataType, typename PrimalElementDataType, typename AdjointElementDataType>
void RunRansEvmKEpsilonElementDataTest(const Variable<TDataType>& rPerturbVariable,
                                       const RansModellingApplicationTestUtilities::ElementDataMethods& rMethodType,
                                       const double Delta,
                                       const double Tolerance)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ElementsContainerType>(
        r_primal_model_part, "Element2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ElementsContainerType>(
        r_adjoint_model_part, "Element2D3N");

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansNutKEpsilonHighReUpdateProcess adjoint_nut_process(adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReUpdateProcess primal_nut_process(primal_model, empty_nut_parameters);

    std::vector<Process*> primal_processes;
    std::vector<Process*> adjoint_processes;

    primal_processes.push_back(&primal_nut_process);
    adjoint_processes.push_back(&adjoint_nut_process);

    auto perturbation_variable =
        RansModellingApplicationTestUtilities::GetPerturbationMethod(rPerturbVariable);

    auto primal_method =
        RansModellingApplicationTestUtilities::GetPrimalMethod<PrimalElementDataType>(rMethodType);
    auto adjoint_method =
        RansModellingApplicationTestUtilities::GetAdjointMethod<TDim, TNumNodes, AdjointElementDataType>(
            rMethodType, rPerturbVariable);

    RansModellingApplicationTestUtilities::RunAdjointElementDataSensitivityTest<TDim, TNumNodes, PrimalElementDataType, AdjointElementDataType>(
        r_primal_model_part, r_adjoint_model_part, primal_processes, adjoint_processes,
        perturbation_variable, primal_method, adjoint_method, Delta, Tolerance);
}

} // namespace RansEvmKEpsilonModel
} // namespace Testing
} // namespace Kratos
#endif
