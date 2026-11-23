---
type: C++ Method
title: RemoveNamespace
resource: sources/code_location.cpp#L85-L90
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/code_location/CodeLocation/ReplaceAll
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/code_location/CodeLocation/CleanFunctionName
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CodeLocation::RemoveNamespace(std::string& FunctionName, const std::string& Namespace)`

# Calls

- [ReplaceAll](../../../../functions/sources/code_location/CodeLocation/ReplaceAll.md)

# Called by

- [CleanFunctionName](../../../../functions/sources/code_location/CodeLocation/CleanFunctionName.md)