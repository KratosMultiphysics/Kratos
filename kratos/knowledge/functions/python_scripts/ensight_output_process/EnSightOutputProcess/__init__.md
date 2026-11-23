---
type: Python Method
title: __init__
resource: python_scripts/ensight_output_process.py#L38-L68
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/EnSightOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/MyPID
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/DeleteDirectoryIfExisting
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __init__(self, model: KratosMultiphysics.Model, settings: KratosMultiphysics.Parameters) -> None:`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [EnSightOutput](../../../../functions/input_output/ensight_output/EnSightOutput/EnSightOutput.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [MyPID](../../../../functions/sources/communicator/Communicator/MyPID.md)
- [DeleteDirectoryIfExisting](../../../../functions/python_scripts/kratos_utilities/DeleteDirectoryIfExisting.md)