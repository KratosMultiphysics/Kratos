---
type: C++ Method
title: GetOutputFileName
resource: input_output/vtk_output.cpp#L224-L274
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/communicator/Communicator/MyPID
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_filesystem/FilesystemExtensions/MPISafeCreateDirectories
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::string VtkOutput::GetOutputFileName(const ModelPart& rModelPart, const bool IsSubModelPart, const std::string& rOutputFilename)`

# Calls

- [MyPID](../../../../functions/sources/communicator/Communicator/MyPID.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [MPISafeCreateDirectories](../../../../functions/sources/kratos_filesystem/FilesystemExtensions/MPISafeCreateDirectories.md)