---
type: C++ Method
title: CloseOutputFile
resource: utilities/timer.cpp#L157-L163
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/timer_process/TimerProcess/ExecuteFinalize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/external_libraries/test_delaunator/validatewithtolerance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/CloseOuputFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int Timer::CloseOutputFile()`

# Called by

- [ExecuteFinalize](../../../../functions/python_scripts/timer_process/TimerProcess/ExecuteFinalize.md)
- [validatewithtolerance](../../../../functions/tests/cpp_tests/external_libraries/test_delaunator/validatewithtolerance.md)
- [CloseOuputFile](../../../../functions/utilities/timer/Timer/CloseOuputFile.md)