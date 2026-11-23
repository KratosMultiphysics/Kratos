---
type: C++ Method
title: UpdateCurrentPosition
resource: utilities/variable_utils.cpp#L273-L286
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python/add_strategies_to_python/MoveMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_variable_utils_to_python/VariableUtilsUpdateCurrentPosition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_variable_utils_to_python/VariableUtilsUpdateCurrentPositionWithVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_variable_utils_to_python/VariableUtilsUpdateCurrentPositionWithVariableAndPosition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildAndSolveLinearizedOnPreviousIteration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_utils/TestVariableUtils/test_UpdateCurrentPosition
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VariableUtils::UpdateCurrentPosition( const ModelPart::NodesContainerType& rNodes, const ArrayVarType& rUpdateVariable, const IndexType BufferPosition )`

# Called by

- [MoveMesh](../../../../functions/python/add_strategies_to_python/MoveMesh.md)
- [VariableUtilsUpdateCurrentPosition](../../../../functions/python/add_variable_utils_to_python/VariableUtilsUpdateCurrentPosition.md)
- [VariableUtilsUpdateCurrentPositionWithVariable](../../../../functions/python/add_variable_utils_to_python/VariableUtilsUpdateCurrentPositionWithVariable.md)
- [VariableUtilsUpdateCurrentPositionWithVariableAndPosition](../../../../functions/python/add_variable_utils_to_python/VariableUtilsUpdateCurrentPositionWithVariableAndPosition.md)
- [BuildAndSolveLinearizedOnPreviousIteration](../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildAndSolveLinearizedOnPreviousIteration.md)
- [test_UpdateCurrentPosition](../../../../functions/tests/test_variable_utils/TestVariableUtils/test_UpdateCurrentPosition.md)