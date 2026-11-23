---
type: Python Method
title: _prepare_test
resource: tests/test_materials_input.py#L18-L39
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_materials_input/TestMaterialsInput/test_input_cpp
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_materials_input/TestMaterialsInput/test_input_with_subproperties_cpp
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_materials_input/TestMaterialsInput/test_input_with_subproperties_cpp_expected_failure
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_materials_input/TestMaterialsInput/test_input_without_tables_and_variables
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _prepare_test(self, input_file = "materials.json"): # Define a Model`

# Calls

- [SetString](../../../../functions/sources/kratos_parameters/Parameters/SetString.md)

# Called by

- [test_input_cpp](../../../../functions/tests/test_materials_input/TestMaterialsInput/test_input_cpp.md)
- [test_input_with_subproperties_cpp](../../../../functions/tests/test_materials_input/TestMaterialsInput/test_input_with_subproperties_cpp.md)
- [test_input_with_subproperties_cpp_expected_failure](../../../../functions/tests/test_materials_input/TestMaterialsInput/test_input_with_subproperties_cpp_expected_failure.md)
- [test_input_without_tables_and_variables](../../../../functions/tests/test_materials_input/TestMaterialsInput/test_input_without_tables_and_variables.md)