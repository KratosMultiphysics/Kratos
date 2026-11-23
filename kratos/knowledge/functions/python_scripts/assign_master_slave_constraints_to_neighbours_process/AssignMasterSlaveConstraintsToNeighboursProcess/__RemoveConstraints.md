---
type: Python Method
title: __RemoveConstraints
resource: python_scripts/assign_master_slave_constraints_to_neighbours_process.py#L138-L141
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/variable_utils/SetFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveMasterSlaveConstraintsFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __RemoveConstraints(self): #Remove master-slave constraints`

# Calls

- [SetFlag](../../../../functions/utilities/variable_utils/SetFlag.md)
- [RemoveMasterSlaveConstraintsFromAllLevels](../../../../functions/sources/model_part/ModelPart/RemoveMasterSlaveConstraintsFromAllLevels.md)