---
type: C++ Function
title: WriteInThisRank
resource: input_output/logger_message.h#L229-L231
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/logger_message/DistributedFilter/WriteFromRank
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/logger_output/LoggerOutput/WriteMessage
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/logger_table_output/LoggerTableOutput/WriteMessage
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool WriteInThisRank() const`

# Calls

- [WriteFromRank](../../../functions/input_output/logger_message/DistributedFilter/WriteFromRank.md)

# Called by

- [WriteMessage](../../../functions/input_output/logger_output/LoggerOutput/WriteMessage.md)
- [WriteMessage](../../../functions/input_output/logger_table_output/LoggerTableOutput/WriteMessage.md)