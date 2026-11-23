---
type: C++ Method
title: WriteMessage
resource: input_output/logger_output.cpp#L53-L92
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/logger_message/WriteInThisRank
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/logger_output/LoggerOutput/SetMessageColor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/logger_message/GetSourceRank
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/logger_output/LoggerOutput/ResetMessageColor
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void LoggerOutput::WriteMessage(LoggerMessage const& TheMessage)`

# Calls

- [WriteInThisRank](../../../../functions/input_output/logger_message/WriteInThisRank.md)
- [SetMessageColor](../../../../functions/input_output/logger_output/LoggerOutput/SetMessageColor.md)
- [GetSourceRank](../../../../functions/input_output/logger_message/GetSourceRank.md)
- [ResetMessageColor](../../../../functions/input_output/logger_output/LoggerOutput/ResetMessageColor.md)