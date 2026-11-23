---
type: C++ Method
title: StructuredMeshGeneratorProcess
resource: processes/structured_mesh_generator_process.cpp#L32-L72
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/ValidateTheDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`StructuredMeshGeneratorProcess::StructuredMeshGeneratorProcess(const GeometryType& rGeometry, ModelPart& rOutputModelPart, Parameters TheParameters) : Process() , mrGeometry(rGeometry) , mrOutputModelPart(rOutputModelPart)`

# Calls

- [ValidateTheDefaultParameters](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/ValidateTheDefaultParameters.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [IsInt](../../../../functions/sources/kratos_parameters/Parameters/IsInt.md)
- [IsVector](../../../../functions/sources/kratos_parameters/Parameters/IsVector.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)