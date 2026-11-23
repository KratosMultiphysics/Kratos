---
type: C++ Method
title: ReduceTemplateArgumentsToFirstN
resource: sources/code_location.cpp#L92-L135
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/code_location/CodeLocation/GetNextPositionSkippingWhiteSpaces
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/code_location/CodeLocation/CleanFunctionName
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CodeLocation::ReduceTemplateArgumentsToFirstN(std::string& FunctionName, const std::string& TemplateName, std::size_t NumberOfArgumentsToKeep)`

# Calls

- [GetNextPositionSkippingWhiteSpaces](../../../../functions/sources/code_location/CodeLocation/GetNextPositionSkippingWhiteSpaces.md)

# Called by

- [CleanFunctionName](../../../../functions/sources/code_location/CodeLocation/CleanFunctionName.md)