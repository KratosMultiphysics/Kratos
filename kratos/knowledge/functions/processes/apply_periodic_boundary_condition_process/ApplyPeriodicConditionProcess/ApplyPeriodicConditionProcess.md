---
type: C++ Method
title: ApplyPeriodicConditionProcess
resource: processes/apply_periodic_boundary_condition_process.cpp#L30-L65
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/RecursivelyValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/RemoveCommonNodesFromSlaveModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/CalculateTransformationMatrix
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/apply_periodic_condition_process/ApplyPeriodicConditionProcessPythonWrapper/__init__
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ApplyPeriodicConditionProcess::ApplyPeriodicConditionProcess(ModelPart &rMasterModelPart, ModelPart &rSlaveModelPart, Parameters Settings) : Process(Flags()), mrMasterModelPart(rMasterModelPart), mrSlaveModelPart(rSlaveModelPart), mParameters(Settings)`

# Calls

- [RecursivelyValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/RecursivelyValidateAndAssignDefaults.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [RemoveCommonNodesFromSlaveModelPart](../../../../functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/RemoveCommonNodesFromSlaveModelPart.md)
- [CalculateTransformationMatrix](../../../../functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/CalculateTransformationMatrix.md)

# Called by

- [__init__](../../../../functions/python_scripts/apply_periodic_condition_process/ApplyPeriodicConditionProcessPythonWrapper/__init__.md)