---
type: C++ Function
title: SplitInternalInterfacesProcess
resource: processes/split_internal_interfaces_process.h#L84-L108
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/HasSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/split_internal_interfaces_process/Factory
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_split_internal_interfaces_process_2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_split_internal_interfaces_process_3D
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`SplitInternalInterfacesProcess( Model& rModel, Parameters rParameters) : Process() , mrModelPart(rModel.GetModelPart(rParameters["model_part_name"].GetString()))`

# Calls

- [GetString](../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [ValidateAndAssignDefaults](../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [HasSubModelPart](../../../functions/sources/model_part/ModelPart/HasSubModelPart.md)

# Called by

- [Factory](../../../functions/python_scripts/split_internal_interfaces_process/Factory.md)
- [test_split_internal_interfaces_process_2D](../../../functions/tests/test_processes/TestProcesses/test_split_internal_interfaces_process_2D.md)
- [test_split_internal_interfaces_process_3D](../../../functions/tests/test_processes/TestProcesses/test_split_internal_interfaces_process_3D.md)