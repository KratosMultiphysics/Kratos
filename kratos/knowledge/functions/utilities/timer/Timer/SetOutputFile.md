---
type: C++ Method
title: SetOutputFile
resource: utilities/timer.cpp#L143-L155
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/timer_process/TimerProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ModelPartIO
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/external_libraries/test_delaunator/validatewithtolerance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/SetOuputFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int Timer::SetOutputFile(std::string const& rOutputFileName)`

# Called by

- [__init__](../../../../functions/python_scripts/timer_process/TimerProcess/__init__.md)
- [ModelPartIO](../../../../functions/sources/model_part_io/ModelPartIO/ModelPartIO.md)
- [validatewithtolerance](../../../../functions/tests/cpp_tests/external_libraries/test_delaunator/validatewithtolerance.md)
- [SetOuputFile](../../../../functions/utilities/timer/Timer/SetOuputFile.md)