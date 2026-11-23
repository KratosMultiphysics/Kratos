---
type: C++ Method
title: CalculateTransformationMatrix
resource: processes/apply_periodic_boundary_condition_process.cpp#L276-L286
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/geometrical_transformation_utilities/CalculateTranslationMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometrical_transformation_utilities/CalculateRotationMatrix
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/ApplyPeriodicConditionProcess
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ApplyPeriodicConditionProcess::CalculateTransformationMatrix()`

# Calls

- [CalculateTranslationMatrix](../../../../functions/utilities/geometrical_transformation_utilities/CalculateTranslationMatrix.md)
- [CalculateRotationMatrix](../../../../functions/utilities/geometrical_transformation_utilities/CalculateRotationMatrix.md)

# Called by

- [ApplyPeriodicConditionProcess](../../../../functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/ApplyPeriodicConditionProcess.md)