---
type: C++ Function
title: printImpl
resource: python/add_logger_to_python.cpp#L46-L97
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python/add_logger_to_python/printDefault
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_logger_to_python/printInfo
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_logger_to_python/printWarning
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_logger_to_python/printDefaultOnAllRanks
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_logger_to_python/printInfoOnAllRanks
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_logger_to_python/printWarningOnAllRanks
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void printImpl(pybind11::args args, pybind11::kwargs kwargs, Logger::Severity severity, bool useKwargLabel, LoggerMessage::DistributedFilter filterOption)`

# Called by

- [printDefault](../../../functions/python/add_logger_to_python/printDefault.md)
- [printInfo](../../../functions/python/add_logger_to_python/printInfo.md)
- [printWarning](../../../functions/python/add_logger_to_python/printWarning.md)
- [printDefaultOnAllRanks](../../../functions/python/add_logger_to_python/printDefaultOnAllRanks.md)
- [printInfoOnAllRanks](../../../functions/python/add_logger_to_python/printInfoOnAllRanks.md)
- [printWarningOnAllRanks](../../../functions/python/add_logger_to_python/printWarningOnAllRanks.md)