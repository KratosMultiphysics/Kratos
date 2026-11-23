---
type: C++ Method
title: EntitiesEraseProcess
resource: processes/entity_erase_process.cpp#L33-L43
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/RecursivelyValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`EntitiesEraseProcess<TEntity>::EntitiesEraseProcess( ModelPart& rModelPart, const Flags Options ) : mrModelPart(rModelPart), mrOptions(Options)`

# Calls

- [RecursivelyValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/RecursivelyValidateAndAssignDefaults.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)