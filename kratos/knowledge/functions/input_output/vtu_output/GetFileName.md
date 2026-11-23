---
type: C++ Function
title: GetFileName
resource: input_output/vtu_output.cpp#L70-L86
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_filesystem/parent_path
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::string GetFileName( const std::string& rBaseName, const std::string& rOutputFileNamePrefix, const bool IsDistributed, const int PartitionIndex)`

# Calls

- [parent_path](../../../functions/sources/kratos_filesystem/parent_path.md)