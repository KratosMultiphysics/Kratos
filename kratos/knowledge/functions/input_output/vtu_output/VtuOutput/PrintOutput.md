---
type: C++ Method
title: PrintOutput
resource: input_output/vtu_output.cpp#L1389-L1393
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/xml_utilities/xml_element/XmlElement/AddAttribute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/MyPID
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/GetEndianness
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/GetBaseName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/TotalProcesses
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_filesystem/parent_path
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtuOutput::PrintOutput(const std::string& rOutputFileNamePrefix)`

# Calls

- [AddAttribute](../../../../functions/utilities/xml_utilities/xml_element/XmlElement/AddAttribute.md)
- [MyPID](../../../../functions/sources/communicator/Communicator/MyPID.md)
- [FullName](../../../../functions/includes/model_part/FullName.md)
- [GetEndianness](../../../../functions/input_output/vtu_output/GetEndianness.md)
- [GetBaseName](../../../../functions/input_output/vtu_output/GetBaseName.md)
- [TotalProcesses](../../../../functions/sources/communicator/Communicator/TotalProcesses.md)
- [parent_path](../../../../functions/sources/kratos_filesystem/parent_path.md)