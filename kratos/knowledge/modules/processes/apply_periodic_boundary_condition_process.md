---
type: C++ Module
title: apply_periodic_boundary_condition_process
resource: processes/apply_periodic_boundary_condition_process.cpp#L1-L304
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/containers-model-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-apply-periodic-boundary-condition-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-binbased-fast-point-locator-conditions-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-geometrical-transformation-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-builtin-timer-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-parallel-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-reduction-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [ApplyPeriodicConditionProcess](../../functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/ApplyPeriodicConditionProcess.md)
- [Create](../../functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/Create.md)
- [RemoveCommonNodesFromSlaveModelPart](../../functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/RemoveCommonNodesFromSlaveModelPart.md)
- [ExecuteInitialize](../../functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/ExecuteInitialize.md)
- [ExecuteInitializeSolutionStep](../../functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/ExecuteInitializeSolutionStep.md)
- [GetDefaultParameters](../../functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/GetDefaultParameters.md)
- [PrintInfo](../../functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/PrintInfo.md)
- [ApplyConstraintsForPeriodicConditions](../../functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/ApplyConstraintsForPeriodicConditions.md)
- [ConstraintSlaveNodeWithConditionForVectorVariable](../../functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/ConstraintSlaveNodeWithConditionForVectorVariable.md)
- [ConstraintSlaveNodeWithConditionForScalarVariable](../../functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/ConstraintSlaveNodeWithConditionForScalarVariable.md)
- [CalculateTransformationMatrix](../../functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/CalculateTransformationMatrix.md)
- [TransformNode](../../functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/TransformNode.md)

# Imports

- `"containers/model.h"`
- `"processes/apply_periodic_boundary_condition_process.h"`
- `"utilities/binbased_fast_point_locator_conditions.h"`
- `"utilities/geometrical_transformation_utilities.h"`
- `"utilities/builtin_timer.h"`
- `"utilities/parallel_utilities.h"`
- `"utilities/reduction_utilities.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)