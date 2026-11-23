---
type: Python Method
title: __init__
resource: python_scripts/timer_process.py#L22-L63
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Timer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/SetPrintIntervalInformation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/SetOutputFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/SetPrintOnScreen
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Start
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __init__(self, Model, settings ):`

# Calls

- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [Timer](../../../../functions/utilities/timer/Timer/Timer.md)
- [SetPrintIntervalInformation](../../../../functions/utilities/timer/Timer/SetPrintIntervalInformation.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [SetOutputFile](../../../../functions/utilities/timer/Timer/SetOutputFile.md)
- [SetPrintOnScreen](../../../../functions/utilities/timer/Timer/SetPrintOnScreen.md)
- [Start](../../../../functions/utilities/timer/Timer/Start.md)