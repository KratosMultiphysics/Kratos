---
type: C++ Method
title: WriteProperties
resource: sources/model_part_io.cpp#L172-L226
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/gid_io/GidIOBase/GetData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_filesystem/remove
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::WriteProperties(PropertiesContainerType const& rThisProperties)`

# Calls

- [GetData](../../../../functions/sources/gid_io/GidIOBase/GetData.md)
- [remove](../../../../functions/sources/kratos_filesystem/remove.md)