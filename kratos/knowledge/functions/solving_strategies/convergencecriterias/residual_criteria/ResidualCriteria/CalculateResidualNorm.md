---
type: C++ Method
title: CalculateResidualNorm
resource: solving_strategies/convergencecriterias/residual_criteria.h#L430-L477
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/convergencecriterias/residual_criteria/ResidualCriteria/IsActiveAndLocalDof
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/convergencecriterias/residual_criteria/ResidualCriteria/PostCriteria
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/convergencecriterias/residual_criteria/ResidualCriteria/InitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void CalculateResidualNorm( ModelPart& rModelPart, TDataType& rResidualSolutionNorm, SizeType& rDofNum, DofsArrayType& rDofSet, const TSystemVectorType& rb )`

# Calls

- [IsActiveAndLocalDof](../../../../../functions/solving_strategies/convergencecriterias/residual_criteria/ResidualCriteria/IsActiveAndLocalDof.md)

# Called by

- [PostCriteria](../../../../../functions/solving_strategies/convergencecriterias/residual_criteria/ResidualCriteria/PostCriteria.md)
- [InitializeSolutionStep](../../../../../functions/solving_strategies/convergencecriterias/residual_criteria/ResidualCriteria/InitializeSolutionStep.md)