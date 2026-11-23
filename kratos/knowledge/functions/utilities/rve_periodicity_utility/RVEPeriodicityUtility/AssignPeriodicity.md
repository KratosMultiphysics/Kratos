---
type: C++ Method
title: AssignPeriodicity
resource: utilities/rve_periodicity_utility.cpp#L23-L45
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void RVEPeriodicityUtility::AssignPeriodicity( ModelPart& rMasterModelPart, ModelPart& rSlaveModelPart, const Matrix& rStrainTensor, const Vector& rDirection, const double SearchTolerance )`

# Calls

- [FullName](../../../../functions/includes/model_part/FullName.md)