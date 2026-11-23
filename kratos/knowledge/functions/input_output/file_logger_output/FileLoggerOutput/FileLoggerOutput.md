---
type: C++ Method
title: FileLoggerOutput
resource: input_output/file_logger_output.cpp#L26-L31
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/logger_output/SetStream
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_file_logger/TestFileLogger/setUpClass
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`FileLoggerOutput::FileLoggerOutput(const std::string& rFileName) : LoggerOutput(), mFileStream(rFileName)`

# Calls

- [SetStream](../../../../functions/input_output/logger_output/SetStream.md)

# Called by

- [setUpClass](../../../../functions/tests/test_file_logger/TestFileLogger/setUpClass.md)