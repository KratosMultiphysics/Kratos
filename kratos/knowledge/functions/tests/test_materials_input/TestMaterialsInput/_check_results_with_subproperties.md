---
type: Python Method
title: _check_results_with_subproperties
resource: tests/test_materials_input.py#L82-L94
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/properties/NumberOfSubproperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/properties/HasSubProperties
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_materials_input/TestMaterialsInput/test_input_with_subproperties_cpp
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _check_results_with_subproperties(self):`

# Calls

- [NumberOfSubproperties](../../../../functions/includes/properties/NumberOfSubproperties.md)
- [HasSubProperties](../../../../functions/includes/properties/HasSubProperties.md)

# Called by

- [test_input_with_subproperties_cpp](../../../../functions/tests/test_materials_input/TestMaterialsInput/test_input_with_subproperties_cpp.md)