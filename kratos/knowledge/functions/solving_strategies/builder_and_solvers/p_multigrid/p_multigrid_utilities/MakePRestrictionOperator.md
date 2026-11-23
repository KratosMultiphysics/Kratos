---
type: C++ Function
title: MakePRestrictionOperator
resource: solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_utilities.hpp#L95-L261
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/Matrix/fill
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/variable_data/IsNotComponent
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/nodal_data/SetSolutionStepData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/GetSolutionStepsData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/SetEquationId
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/FixDof
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MakePRestrictionOperator(const Geometry<TNode>& rGeometry, TOutputIterator itOutput)`

# Calls

- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [fill](../../../../../functions/includes/amatrix_interface/Matrix/fill.md)
- [GetVariable](../../../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)
- [IsNotComponent](../../../../../functions/containers/variable_data/IsNotComponent.md)
- [SetSolutionStepData](../../../../../functions/containers/nodal_data/SetSolutionStepData.md)
- [GetSolutionStepsData](../../../../../functions/includes/dof/Dof/GetSolutionStepsData.md)
- [SetEquationId](../../../../../functions/includes/dof/Dof/SetEquationId.md)
- [FixDof](../../../../../functions/includes/dof/Dof/FixDof.md)