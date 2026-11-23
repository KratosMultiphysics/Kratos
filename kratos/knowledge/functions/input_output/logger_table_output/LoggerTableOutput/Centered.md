---
type: C++ Method
title: Centered
resource: input_output/logger_table_output.cpp#L212-L217
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/logger_table_output/LoggerTableOutput/WriteHeader
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/logger_table_output/LoggerTableOutput/WriteTableLine
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::string LoggerTableOutput::Centered(int width, const std::string& str)`

# Called by

- [WriteHeader](../../../../functions/input_output/logger_table_output/LoggerTableOutput/WriteHeader.md)
- [WriteTableLine](../../../../functions/input_output/logger_table_output/LoggerTableOutput/WriteTableLine.md)