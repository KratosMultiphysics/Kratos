---
type: Python Method
title: __init__
resource: python_scripts/vtk_output_process.py#L14-L38
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddEmptyValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/RemoveValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/VtkOutput
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
- [AddEmptyValue](../../../../functions/sources/kratos_parameters/Parameters/AddEmptyValue.md)
- [SetBool](../../../../functions/sources/kratos_parameters/Parameters/SetBool.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [RemoveValue](../../../../functions/sources/kratos_parameters/Parameters/RemoveValue.md)
- [VtkOutput](../../../../functions/input_output/vtk_output/VtkOutput/VtkOutput.md)
- [MyPID](../../../../functions/sources/communicator/Communicator/MyPID.md)
- [DeleteDirectoryIfExisting](../../../../functions/python_scripts/kratos_utilities/DeleteDirectoryIfExisting.md)