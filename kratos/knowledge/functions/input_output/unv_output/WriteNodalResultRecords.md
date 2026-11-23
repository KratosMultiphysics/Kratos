---
type: C++ Function
title: WriteNodalResultRecords
resource: input_output/unv_output.h#L225-L298
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/unv_output/as_integer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/unv_output/GetUnvVariableName
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/unv_output/UnvOutput/WriteNodalResults
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void WriteNodalResultRecords(const TVariablebleType& rVariable, const int numComponents, const double timeStep)`

# Calls

- [as_integer](../../../functions/input_output/unv_output/as_integer.md)
- [GetUnvVariableName](../../../functions/input_output/unv_output/GetUnvVariableName.md)

# Called by

- [WriteNodalResults](../../../functions/input_output/unv_output/UnvOutput/WriteNodalResults.md)