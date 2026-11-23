---
type: C++ Method
title: Execute
resource: processes/check_same_modelpart_using_skin_distance_process.cpp#L31-L151
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/project/Project/GetModel
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/model/Model/DeleteModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CheckSameModelPartUsingSkinDistanceProcess<TDim>::Execute()`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [GetModel](../../../../functions/python_scripts/project/Project/GetModel.md)
- [SetInt](../../../../functions/sources/kratos_parameters/Parameters/SetInt.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [SetString](../../../../functions/sources/kratos_parameters/Parameters/SetString.md)
- [AuxiliarModelPartUtilities](../../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)
- [DeleteModelPart](../../../../functions/containers/model/Model/DeleteModelPart.md)