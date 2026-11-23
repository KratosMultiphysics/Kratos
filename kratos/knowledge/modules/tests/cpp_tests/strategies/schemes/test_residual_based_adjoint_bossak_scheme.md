---
type: C++ Module
title: test_residual_based_adjoint_bossak_scheme
resource: tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme.cpp#L1-L1104
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/algorithm
    resolved_by: tree-sitter
    confidence: exact
  - target: external/exception
    resolved_by: tree-sitter
    confidence: exact
  - target: external/sstream
    resolved_by: tree-sitter
    confidence: exact
  - target: external/testing-testing-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-model-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-define-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-model-part-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/spaces-ublas-space-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-skyline-lu-custom-scalar-solver-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-schemes-residual-based-bossak-displacement-scheme-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-schemes-residual-based-adjoint-bossak-scheme-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-strategies-residualbased-newton-raphson-strategy-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-strategies-residualbased-linear-strategy-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-convergencecriterias-residual-criteria-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-indirect-scalar-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/response-functions-adjoint-response-function-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-sensitivity-builder-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-adjoint-extensions-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [PrimalResults](../../../../../classes/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalResults.md)
- [PrimalStrategy](../../../../../classes/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalStrategy.md)
- [PrimalStrategy](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalStrategy/PrimalStrategy.md)
- [Initialize](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalStrategy/Initialize.md)
- [Solve](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalStrategy/Solve.md)
- [AdjointStrategy](../../../../../classes/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointStrategy.md)
- [AdjointStrategy](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointStrategy/AdjointStrategy.md)
- [Initialize](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointStrategy/Initialize.md)
- [Solve](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointStrategy/Solve.md)
- [PrimalElement](../../../../../classes/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalElement.md)
- [CreateEntity](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalElement/CreateEntity.md)
- [PrimalElement](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalElement/PrimalElement.md)
- [EquationIdVector](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalElement/EquationIdVector.md)
- [GetDofList](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalElement/GetDofList.md)
- [GetValuesVector](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalElement/GetValuesVector.md)
- [GetFirstDerivativesVector](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalElement/GetFirstDerivativesVector.md)
- [GetSecondDerivativesVector](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalElement/GetSecondDerivativesVector.md)
- [CalculateLocalSystem](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalElement/CalculateLocalSystem.md)
- [CalculateLeftHandSide](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalElement/CalculateLeftHandSide.md)
- [CalculateRightHandSide](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalElement/CalculateRightHandSide.md)
- [CalculateMassMatrix](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalElement/CalculateMassMatrix.md)
- [CalculateDampingMatrix](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalElement/CalculateDampingMatrix.md)
- [PrimalCondition](../../../../../classes/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalCondition.md)
- [CreateEntity](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalCondition/CreateEntity.md)
- [PrimalCondition](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalCondition/PrimalCondition.md)
- [EquationIdVector](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalCondition/EquationIdVector.md)
- [GetDofList](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalCondition/GetDofList.md)
- [GetValuesVector](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalCondition/GetValuesVector.md)
- [GetFirstDerivativesVector](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalCondition/GetFirstDerivativesVector.md)
- [GetSecondDerivativesVector](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalCondition/GetSecondDerivativesVector.md)
- [CalculateLocalSystem](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalCondition/CalculateLocalSystem.md)
- [CalculateLeftHandSide](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalCondition/CalculateLeftHandSide.md)
- [CalculateRightHandSide](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalCondition/CalculateRightHandSide.md)
- [CalculateMassMatrix](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalCondition/CalculateMassMatrix.md)
- [CalculateDampingMatrix](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalCondition/CalculateDampingMatrix.md)
- [AdjointElement](../../../../../classes/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointElement.md)
- [ThisExtensions](../../../../../classes/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ThisExtensions.md)
- [ThisExtensions](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ThisExtensions/ThisExtensions.md)
- [GetFirstDerivativesVector](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ThisExtensions/GetFirstDerivativesVector.md)
- [GetSecondDerivativesVector](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ThisExtensions/GetSecondDerivativesVector.md)
- [GetAuxiliaryVector](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ThisExtensions/GetAuxiliaryVector.md)
- [GetFirstDerivativesVariables](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ThisExtensions/GetFirstDerivativesVariables.md)
- [GetSecondDerivativesVariables](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ThisExtensions/GetSecondDerivativesVariables.md)
- [GetAuxiliaryVariables](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ThisExtensions/GetAuxiliaryVariables.md)
- [CreateEntity](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointElement/CreateEntity.md)
- [AdjointElement](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointElement/AdjointElement.md)
- [EquationIdVector](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointElement/EquationIdVector.md)
- [GetDofList](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointElement/GetDofList.md)
- [GetValuesVector](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointElement/GetValuesVector.md)
- [GetFirstDerivativesVector](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointElement/GetFirstDerivativesVector.md)
- [GetSecondDerivativesVector](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointElement/GetSecondDerivativesVector.md)
- [CalculateLeftHandSide](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointElement/CalculateLeftHandSide.md)
- [CalculateFirstDerivativesLHS](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointElement/CalculateFirstDerivativesLHS.md)
- [CalculateSecondDerivativesLHS](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointElement/CalculateSecondDerivativesLHS.md)
- [CalculateSensitivityMatrix](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointElement/CalculateSensitivityMatrix.md)
- [AdjointCondition](../../../../../classes/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointCondition.md)
- [ThisExtensions](../../../../../classes/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ThisExtensions-2.md)
- [ThisExtensions](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ThisExtensions/ThisExtensions-2.md)
- [GetFirstDerivativesVector](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ThisExtensions/GetFirstDerivativesVector-2.md)
- [GetSecondDerivativesVector](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ThisExtensions/GetSecondDerivativesVector-2.md)
- [GetAuxiliaryVector](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ThisExtensions/GetAuxiliaryVector-2.md)
- [GetFirstDerivativesVariables](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ThisExtensions/GetFirstDerivativesVariables-2.md)
- [GetSecondDerivativesVariables](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ThisExtensions/GetSecondDerivativesVariables-2.md)
- [GetAuxiliaryVariables](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ThisExtensions/GetAuxiliaryVariables-2.md)
- [CreateEntity](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointCondition/CreateEntity.md)
- [AdjointCondition](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointCondition/AdjointCondition.md)
- [EquationIdVector](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointCondition/EquationIdVector.md)
- [GetDofList](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointCondition/GetDofList.md)
- [GetValuesVector](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointCondition/GetValuesVector.md)
- [GetFirstDerivativesVector](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointCondition/GetFirstDerivativesVector.md)
- [GetSecondDerivativesVector](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointCondition/GetSecondDerivativesVector.md)
- [CalculateLeftHandSide](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointCondition/CalculateLeftHandSide.md)
- [CalculateFirstDerivativesLHS](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointCondition/CalculateFirstDerivativesLHS.md)
- [CalculateSecondDerivativesLHS](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointCondition/CalculateSecondDerivativesLHS.md)
- [CalculateSensitivityMatrix](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointCondition/CalculateSensitivityMatrix.md)
- [ResponseFunction](../../../../../classes/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ResponseFunction.md)
- [ResponseFunction](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ResponseFunction/ResponseFunction.md)
- [CalculateGradient](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ResponseFunction/CalculateGradient.md)
- [CalculateFirstDerivativesGradient](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ResponseFunction/CalculateFirstDerivativesGradient.md)
- [CalculateSecondDerivativesGradient](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ResponseFunction/CalculateSecondDerivativesGradient.md)
- [CalculatePartialSensitivity](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ResponseFunction/CalculatePartialSensitivity.md)
- [CalculatePartialSensitivity](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ResponseFunction/CalculatePartialSensitivity-2.md)
- [CalculateGradient](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ResponseFunction/CalculateGradient-2.md)
- [CalculateFirstDerivativesGradient](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ResponseFunction/CalculateFirstDerivativesGradient-2.md)
- [CalculateSecondDerivativesGradient](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ResponseFunction/CalculateSecondDerivativesGradient-2.md)
- [CalculatePartialSensitivity](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ResponseFunction/CalculatePartialSensitivity-3.md)
- [CalculatePartialSensitivity](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ResponseFunction/CalculatePartialSensitivity-4.md)
- [CalculateValue](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ResponseFunction/CalculateValue.md)
- [PrimalResults](../../../../../classes/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalResults-2.md)
- [StoreCurrentSolutionStep](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalResults/StoreCurrentSolutionStep.md)
- [LoadCurrentSolutionStep](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalResults/LoadCurrentSolutionStep.md)
- [InitializePrimalModelPart](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/InitializePrimalModelPart.md)
- [InitializeAdjointModelPart](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/InitializeAdjointModelPart.md)
- [RunAdjointSensitivityTest](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/RunAdjointSensitivityTest.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/KRATOS_TEST_CASE_IN_SUITE-2.md)

# Imports

- `<algorithm>`
- `<exception>`
- `<sstream>`
- `"testing/testing.h"`
- `"containers/model.h"`
- `"includes/define.h"`
- `"includes/model_part.h"`
- `"spaces/ublas_space.h"`
- `"linear_solvers/skyline_lu_custom_scalar_solver.h"`
- `"solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp"`
- `"solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h"`
- `"solving_strategies/strategies/residualbased_newton_raphson_strategy.h"`
- `"solving_strategies/strategies/residualbased_linear_strategy.h"`
- `"solving_strategies/convergencecriterias/residual_criteria.h"`
- `"utilities/indirect_scalar.h"`
- `"response_functions/adjoint_response_function.h"`
- `"utilities/sensitivity_builder.h"`
- `"utilities/adjoint_extensions.h"`

# Member of

- [KratosMultiphysics](../../../../../packages/KratosMultiphysics.md)