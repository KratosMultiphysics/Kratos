---
type: Python Method
title: _auxiliary_test_function_ExplicitBuilder
resource: tests/test_factories.py#L6-L9
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_factories/TestFactories/test_ExplicitBuilder
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _auxiliary_test_function_ExplicitBuilder(self, settings, name):`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)

# Called by

- [test_ExplicitBuilder](../../../../functions/tests/test_factories/TestFactories/test_ExplicitBuilder.md)