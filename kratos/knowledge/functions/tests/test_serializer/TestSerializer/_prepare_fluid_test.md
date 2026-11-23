---
type: Python Method
title: _prepare_fluid_test
resource: tests/test_serializer.py#L17-L41
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_serializer/TestSerializer/test_serializer_fluid_analysis
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _prepare_fluid_test(self): # Define a model and load the parameters`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [SetString](../../../../functions/sources/kratos_parameters/Parameters/SetString.md)

# Called by

- [test_serializer_fluid_analysis](../../../../functions/tests/test_serializer/TestSerializer/test_serializer_fluid_analysis.md)