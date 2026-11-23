---
type: C++ Method
title: ValidateTheDefaultParameters
resource: processes/structured_mesh_generator_process.cpp#L127-L143
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/IsInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/StructuredMeshGeneratorProcess
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void StructuredMeshGeneratorProcess::ValidateTheDefaultParameters(Parameters TheParameters)`

# Calls

- [IsInt](../../../../functions/sources/kratos_parameters/Parameters/IsInt.md)
- [AddInt](../../../../functions/sources/kratos_parameters/Parameters/AddInt.md)
- [IsVector](../../../../functions/sources/kratos_parameters/Parameters/IsVector.md)
- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)

# Called by

- [StructuredMeshGeneratorProcess](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/StructuredMeshGeneratorProcess.md)