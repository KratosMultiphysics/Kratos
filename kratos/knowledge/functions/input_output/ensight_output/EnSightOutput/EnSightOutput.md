---
type: C++ Method
title: EnSightOutput
resource: input_output/ensight_output.cpp#L68-L158
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GlobalNumberOfElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GlobalNumberOfConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/UpdatePartData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/WriteGeometryFile
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/ensight_output_process/EnSightOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`EnSightOutput::EnSightOutput(ModelPart& rModelPart, Parameters ThisParameters) : mrModelPart(rModelPart), mOutputSettings(ThisParameters)`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GlobalNumberOfElements](../../../../functions/sources/communicator/Communicator/GlobalNumberOfElements.md)
- [GlobalNumberOfConditions](../../../../functions/sources/communicator/Communicator/GlobalNumberOfConditions.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [UpdatePartData](../../../../functions/input_output/ensight_output/EnSightOutput/UpdatePartData.md)
- [WriteGeometryFile](../../../../functions/input_output/ensight_output/EnSightOutput/WriteGeometryFile.md)

# Called by

- [__init__](../../../../functions/python_scripts/ensight_output_process/EnSightOutputProcess/__init__.md)