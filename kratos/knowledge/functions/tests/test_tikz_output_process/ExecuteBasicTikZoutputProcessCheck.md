---
type: Python Function
title: ExecuteBasicTikZoutputProcessCheck
resource: tests/test_tikz_output_process.py#L86-L127
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tikz_output_process/SetupTikZOutputProcess
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_tikz_output_process/TestTikZOutputProcess/test_ascii_tikz_output_2D
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ExecuteBasicTikZoutputProcessCheck():`

# Calls

- [SetString](../../../functions/sources/kratos_parameters/Parameters/SetString.md)
- [SetupTikZOutputProcess](../../../functions/tests/test_tikz_output_process/SetupTikZOutputProcess.md)

# Called by

- [test_ascii_tikz_output_2D](../../../functions/tests/test_tikz_output_process/TestTikZOutputProcess/test_ascii_tikz_output_2D.md)