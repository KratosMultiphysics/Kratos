---
type: C++ Method
title: StripLabels
resource: input_output/logger_table_output.cpp#L218-L225
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/logger_table_output/LoggerTableOutput/WriteMessage
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void LoggerTableOutput::StripLabels(std::string& inLabel, std::string& outGlobalLabel, std::string& outColumnLabel)`

# Called by

- [WriteMessage](../../../../functions/input_output/logger_table_output/LoggerTableOutput/WriteMessage.md)