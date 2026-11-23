---
type: C++ Method
title: Execute
resource: processes/entity_erase_process.cpp#L114-L129
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/variable_utils/SetFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveConditionsFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveConditions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void EntitiesEraseProcess<Condition>::Execute()`

# Calls

- [SetFlag](../../../../functions/utilities/variable_utils/SetFlag.md)
- [RemoveConditionsFromAllLevels](../../../../functions/sources/model_part/ModelPart/RemoveConditionsFromAllLevels.md)
- [RemoveConditions](../../../../functions/sources/model_part/ModelPart/RemoveConditions.md)