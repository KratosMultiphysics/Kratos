---
type: C++ Method
title: WriteMessage
resource: input_output/logger_table_output.cpp#L133-L159
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/logger_table_output/LoggerTableOutput/StripLabels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/logger_message/WriteInThisRank
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/logger_table_output/LoggerTableOutput/WriteTableLine
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void LoggerTableOutput::WriteMessage(LoggerMessage const& TheMessage)`

# Calls

- [StripLabels](../../../../functions/input_output/logger_table_output/LoggerTableOutput/StripLabels.md)
- [WriteInThisRank](../../../../functions/input_output/logger_message/WriteInThisRank.md)
- [WriteTableLine](../../../../functions/input_output/logger_table_output/LoggerTableOutput/WriteTableLine.md)