---
type: C++ Method
title: Write
resource: utilities/profiler.h#L290-L354
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/profiler/Profiler/Aggregate
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddEmptyArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/code_location/CodeLocation/GetLineNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/code_location/CodeLocation/CleanFunctionName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Write(std::ostream& rStream) const`

# Calls

- [Aggregate](../../../../functions/utilities/profiler/Profiler/Aggregate.md)
- [AddString](../../../../functions/sources/kratos_parameters/Parameters/AddString.md)
- [AddInt](../../../../functions/sources/kratos_parameters/Parameters/AddInt.md)
- [AddValue](../../../../functions/sources/kratos_parameters/Parameters/AddValue.md)
- [AddEmptyArray](../../../../functions/sources/kratos_parameters/Parameters/AddEmptyArray.md)
- [GetLineNumber](../../../../functions/sources/code_location/CodeLocation/GetLineNumber.md)
- [CleanFunctionName](../../../../functions/sources/code_location/CodeLocation/CleanFunctionName.md)
- [PrettyPrintJsonString](../../../../functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString.md)