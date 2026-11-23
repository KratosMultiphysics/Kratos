---
type: C++ Method
title: WriteTableLine
resource: input_output/logger_table_output.cpp#L193-L211
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/logger_table_output/LoggerTableOutput/Centered
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/logger_table_output/LoggerTableOutput/WriteMessage
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void LoggerTableOutput::WriteTableLine()`

# Calls

- [Centered](../../../../functions/input_output/logger_table_output/LoggerTableOutput/Centered.md)

# Called by

- [WriteMessage](../../../../functions/input_output/logger_table_output/LoggerTableOutput/WriteMessage.md)