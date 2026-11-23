---
type: C++ Method
title: LoggerTableOutput
resource: input_output/logger_table_output.cpp#L27-L52
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`LoggerTableOutput::LoggerTableOutput(Parameters rSettings): LoggerOutput(mMyFileStream), mCurrentColumnIndex(0), mColumnsHeaders({}), mColumnsWidth({}), mColumnsTexts({}), mColumnsLabels({}), mFileHeader("")`

# Calls

- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)