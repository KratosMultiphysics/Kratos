---
type: C++ Method
title: GetOutputFileName
resource: input_output/ensight_output.cpp#L2142-L2193
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

`std::string EnSightOutput::GetOutputFileName( const std::string& rFileLabel, const std::string& rFileExtension ) const`

# Calls

- [MyPID](../../../../functions/sources/communicator/Communicator/MyPID.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [MPISafeCreateDirectories](../../../../functions/sources/kratos_filesystem/FilesystemExtensions/MPISafeCreateDirectories.md)